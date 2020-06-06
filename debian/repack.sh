#! /bin/bash
#
# debian/repack
# Part of the Debian package ‘libjs-zxcvbn’.
#
# Copyright © 2013–2014 Ben Finney <ben+debian@benfinney.id.au>
# This is free software; see the end of this file for license terms.

# Convert the pristine upstream source to the Debian upstream source.
#
# This program is designed for use with the ‘uscan(1)’ tool, as the
# “action” parameter for the ‘debian/watch’ configuration file.

set -o errexit
set -o errtrace
set -o pipefail
set -o nounset

program_dir="$(dirname "$(realpath --strip "$0")")"
printf "program_dir: ${program_dir}\n"

source "${program_dir}"/source_package_build.bash

function usage() {
	local progname=$(basename $0)
	printf "$progname --upstream-version VERSION FILENAME\n"
}

if [ $# -ne 3 ] ; then
	usage
	printf "Error with the number of arguments: $# instead of 3\n"
	exit 1
fi

upstream_version="$2"
downloaded_file="$3"

target_filename="${upstream_tarball_basename}.tar.gz"
printf "target_filename: ${target_filename}\n"

target_working_file="${working_dir}/${target_filename}"
target_file="$(dirname "${downloaded_file}")/${target_filename}"
printf "target_file: ${target_file}\n"

repack_dir="${working_dir}/${upstream_dirname}"
printf "repack_dir: ${repack_dir}\n"

printf "Unpacking pristine upstream source ‘${downloaded_file}’:\n"

extract_tarball_to_working_dir "${downloaded_file}"

upstream_source_dirname=$(ls -1 "${working_dir}")
upstream_source_dir="${working_dir}/${upstream_source_dirname}"

printf "Repackaging upstream source from ‘${upstream_source_dir}’ to ‘${repack_dir}’:\n"

mv "${upstream_source_dir}" "${repack_dir}"

printf "Removing non-DFSG-free files:\n"

undesirable_fileglobs=(
	IsoSpecPy/IsoSpecPy/prebuilt-libIsoSpec++2.0.1-x32.dll
	IsoSpecPy/IsoSpecPy/prebuilt-libIsoSpec++2.0.1-x64.dll
	man/html
	man/man
	man/latex
	)

for fileglob in "${undesirable_fileglobs[@]}" ; do
	rm -rfv "${repack_dir}"/$fileglob
done

printf "Rebuilding DFSG-free upstream source tarball:\n"

archive_working_dirname_to_tarball "${upstream_dirname}" "${target_working_file}"

printf "Moving completed upstream tarball to ‘${target_file}’:\n"

rm -v "${downloaded_file}"
mv "${target_working_file}" "${target_file}"

printf "Done.\n"


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
# indent-tabs-mode: nil
# End:
# vim: fileencoding=utf-8 filetype=bash :

