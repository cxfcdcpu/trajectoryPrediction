#!/bin/bash
var=0
while [ 1 ]; do
  ./client
  var=$((var+1))
  sleep 3
done
