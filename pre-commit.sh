#!/bin/sh
#
nosetests

rc=$?
if [ "$rc" != "0" ]; then
  echo -n "Not all tests passed. Continue commit? (y/n):"
  read response
  if [[ "$response" == "n" ]]; then
    exit 1
  fi
fi

git diff --cached --name-status | while read st file; do
  if [ "$st" == 'D' ]; then continue; fi
  echo "Checking $file with pylint"
  if [[ "$file" !=  ".py$" && ! $(pylint -E --rcfile pylintrc "$file") ]]; then
    echo "pylint check failed for file $file"
    exit 1
  fi
done

find . -type f -name "*.py" -print0 | xargs -0 pylint --rcfile=pylintrc | grep "Your code has been rated"

