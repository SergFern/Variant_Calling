#!/bin/bash

version="program version #0.1
Written by bioinfo-serg."

usage="Usage: $0 [OPTION]... [FILE]...
Description

Mandatory arguments to long options are mandatory for short options too.

  -short_option, --long_option	description
      --help        display this help and exit
      --version     display version information and exit

With no FILE read standard input.

Report bugs to <your therapist>."

case $1 in
--help)    exec echo "$usage";;
--version) exec echo "$version";;
esac
