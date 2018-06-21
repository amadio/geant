#!/bin/bash
for file in mt*.sh
do
  bash FullCMS.run $file
done
