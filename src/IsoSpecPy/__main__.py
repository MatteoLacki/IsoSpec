from . import __version__

if __name__ == "__main__":
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description="IsoSpecPy: Python interface to IsoSpec++ library, a fast and memory-efficient library for computing isotopic distributions.")
    parser.add_argument('--version', '-v', action='version', version=__version__)
    parser.add_argument('--libpath', action='store_true',
                        help='Print the path to the loaded C++ library and exit.')
    parser.add_argument('--include', action='store_true',
                        help='Print the include path for the headers of the C++ library and exit.')
    args = parser.parse_args()

    if args.libpath:
        try:
            from .isoFFI import IsoFFI
            ffi = IsoFFI()
            print(Path(ffi.libpath).resolve())
        except ImportError as e:
            print(f"Error loading IsoSpecPy: {e}")
        exit(0)

    if args.include:
        print(Path(__file__).parent.resolve())
