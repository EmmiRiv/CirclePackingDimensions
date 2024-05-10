#!/bin/bash

for ((i = 20 ; i < 105 ; i+=5 )); do julia bai_finch.jl $i; done
