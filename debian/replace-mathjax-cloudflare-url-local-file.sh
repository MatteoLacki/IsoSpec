#!/bin/sh

cd man/html

for file in *.html
do
	grep -i "src=\"https://cdnjs.cloudflare.com/ajax/libs/mathjax/[0-9].[0-9].[0-9]/MathJax.js\"" ${file}

	if [ "$?" = "0" ]
	then
		sed -i 's|src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/\([0-9].[0-9].[0-9]\)/MathJax.js\"|src="file:///usr/share/javascript/mathjax/MathJax.js"|' ${file} 
	fi
done
