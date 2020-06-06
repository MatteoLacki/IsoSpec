# debian/source_package_build.bash
# Part of the Debian package ‘libjs-zxcvbn’.
#
# Copyright © 2013–2014 Ben Finney <ben+debian@benfinney.id.au>
# This is free software; see the end of this file for license terms.

# Common code for building Debian upstream source package.

working_dir="$(mktemp -d -t)"

exit_sigspecs="ERR EXIT SIGTERM SIGHUP SIGINT SIGQUIT"

function cleanup_exit() {
    exit_status=$?
    trap - $exit_sigspecs

    rm -rf "${working_dir}"
    printf "Cleaned up working directory ‘${working_dir}’\n"

    exit $exit_status
}
trap cleanup_exit $exit_sigspecs

package_name=$(dpkg-parsechangelog --show-field 'Source')
release_version=$(dpkg-parsechangelog --show-field 'Version')
upstream_version=$(printf "${release_version}" \
        | sed -e 's/^[[:digit:]]\+://' -e 's/[-][^-]\+$//')
upstream_dirname="${package_name}-${upstream_version}.orig"
upstream_tarball_basename="${package_name}_${upstream_version}.orig"

function extract_tarball_to_working_dir() {
    # Extract the specified tarball to the program's working directory.
    local tarball="$1"
    tar -xzf "${tarball}" --directory "${working_dir}"
}

function archive_working_dirname_to_tarball() {
    # Archive the specified directory, relative to the working directory,
    # to a new tarball of the specified name.
    local source_dirname="$1"
    local tarball="$2"
    GZIP="--best" tar --directory "${working_dir}" -czf "${tarball}" "${source_dirname}"
}


# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# “Software”), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# The Software is provided “as is”, without warranty of any kind,
# express or implied, including but not limited to the warranties of
# merchantability, fitness for a particular purpose and noninfringement.
# In no event shall the authors or copyright holders be liable for any
# claim, damages or other liability, whether in an action of contract,
# tort or otherwise, arising from, out of or in connection with the
# Software or the use or other dealings in the Software.


# Local variables:
# coding: utf-8
# mode: sh
# End:
# vim: fileencoding=utf-8 filetype=bash :
