// Builds the IsoSpec C++ core from src/IsoSpec++ together with the cxx bridge,
// the way src/IsoSpec++/IsaDispatch.cmake does for CMake builds:
//
// - unity-build.cpp (every library source, isa_dispatch.cpp included) and the
//   bridge shim, compiled with the target's baseline code generation;
// - the four ISA kernel objects, isa_kernels_{baseline,avx,v3,v4}.cpp, each
//   with its own instruction-set flags -- when the compiler accepts them. A
//   kernel object compiled without its flags publishes an empty table (the
//   sources check their own macros), so on aarch64, with MSVC, or with a
//   compiler lacking a level, that level is simply absent at run time.
//
// All objects share every define: the kernels instantiate the same templates
// as the main build, and a mismatch would be an ODR violation.
//
// The C++ sources come from src/IsoSpec++ in the repository and in the
// published crate alike: Cargo.toml's `include` list ships them.

use std::env;
use std::path::PathBuf;

const CXX_DIR: &str = "src/IsoSpec++";

fn main() {
    let root = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let cxx_dir = root.join(CXX_DIR);

    println!("cargo:rerun-if-changed=rust/src/ffi.rs");
    println!("cargo:rerun-if-changed=rust/shim");
    println!("cargo:rerun-if-changed={CXX_DIR}");
    println!("cargo:rerun-if-env-changed=DOCS_RS");
    println!("cargo:rerun-if-changed=VERSION");

    let version = std::fs::read_to_string(root.join("VERSION")).expect("reading VERSION");
    println!("cargo:rustc-env=ISOSPEC_LIBRARY_VERSION={}", version.trim());

    // docs.rs: rustdoc never links, and its compiler may lack C++20's
    // <experimental/simd> support; documentation needs none of the C++.
    if env::var_os("DOCS_RS").is_some() {
        return;
    }

    let base = configure(cxx_build::bridge("rust/src/ffi.rs"), &cxx_dir);
    let is_msvc = base.get_compiler().is_like_msvc();
    let is_x86 = matches!(
        env::var("CARGO_CFG_TARGET_ARCH").as_deref(),
        Ok("x86_64") | Ok("x86")
    );

    // The kernels: compiled per level with that level's flags, then archived
    // together with everything else, so that one static library resolves all
    // of isa_dispatch's references whatever order the linker reads things in.
    let levels: [(&str, &[&str]); 4] = [
        ("baseline", &[]),
        ("avx", &["-mavx"]),
        ("v3", &["-mavx2", "-mfma"]),
        (
            "v4",
            &["-mavx512f", "-mavx512bw", "-mavx512dq", "-mavx512vl"],
        ),
    ];
    let mut kernel_objects = Vec::new();
    for (level, flags) in levels {
        let mut build = configure(cc::Build::new(), &cxx_dir);
        // MSVC has no <experimental/simd>: its kernels are scalar at every
        // level, so the flags would buy nothing. Off x86 the flags do not exist.
        let wanted = !is_msvc && is_x86 && !flags.is_empty();
        if wanted
            && flags
                .iter()
                .all(|f| build.is_flag_supported(f).unwrap_or(false))
        {
            for f in flags {
                build.flag(f);
            }
        }
        kernel_objects.extend(
            build
                .file(cxx_dir.join(format!("isa_kernels_{level}.cpp")))
                .compile_intermediates(),
        );
    }

    let mut base = base;
    base.file(cxx_dir.join("unity-build.cpp"))
        .file("rust/shim/shim.cc")
        .objects(kernel_objects)
        .compile("isospec_cxx");
}

fn configure(mut build: cc::Build, cxx_dir: &std::path::Path) -> cc::Build {
    build
        .cpp(true)
        .std("c++20")
        .include(cxx_dir)
        .include("rust/shim")
        .define("ISOSPEC_ISA_DISPATCH", "1")
        // The library's own warnings are the C++ build's business; a crate
        // user cannot act on them.
        .warnings(false);
    // The C++ asserts follow Cargo's debug assertions: off in release, as in
    // the library's own release builds (cc does not define NDEBUG itself, and
    // the asserts sit in the per-configuration loops).
    if env::var_os("CARGO_CFG_DEBUG_ASSERTIONS").is_none() {
        build.define("NDEBUG", None);
    }
    if build.get_compiler().is_like_msvc() {
        build.flag("/utf-8").flag("/bigobj").flag("/Zc:__cplusplus");
    }
    build
}
