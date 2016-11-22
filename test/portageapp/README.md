# PortageApp Testing

## How Tests Work

CTest runs the bash `*.sh` files, which run portage app in the various modes
that we want to test (cell/node remap, serial or parallel) and saves the final
field to a text file `fieldn.txt` where `n` is the number of the example in
the portageapp. This `fieldn.txt` is then compared with `field_goldn.txt` using
the comparator app `apptest_cmp`.

## How to Generate Gold Standard Files

Add a new test, without the gold file at first. Run it, it will fail, because
the gold file is not there. Copy the `fieldn.txt` from the build directory into
this test directory and rename it to `field_goldn.txt`.
