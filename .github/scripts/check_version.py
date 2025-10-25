import subprocess
import sys
import tomllib

def git_version():
    version = subprocess.check_output(
        ["git", "describe", "--tags", "--abbrev=0"],
        stderr=subprocess.STDOUT
    ).strip().decode('utf-8')
    return version

def pyproject_version():
    with open("pyproject.toml", "rb") as f:
        pyproject_data = tomllib.load(f)
    version = pyproject_data["project"]["version"]
    return version

if __name__ == "__main__":
    assert git_version() == "v"+pyproject_version(), \
        f"Version mismatch: git version '{git_version()}' != pyproject.toml version '{pyproject_version()}'"
    print("Version check passed:", git_version())