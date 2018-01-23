/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

# PortageApp Testing

## How Tests Work

CTest runs the bash `*.sh` files, which run portage app in the various
modes that we want to test (cell/node remap, serial or parallel) and
saves the final field to a text file `some_name.txt` where
`some_name.txt` is specified using the --results_file argument to the
app. This file is then compared with `GOLD_some_name.txt` using the
comparator app `apptest_cmp`.

## How to Generate Gold Standard Files

Add a new test, with an empty gold file at first. Run it, it will fail, because
the gold file is empty. Copy the `some_name.txt` from the build directory into
this test directory and rename it to `GOLD_some_name.txt`.
