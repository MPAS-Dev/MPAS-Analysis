# License Header Tools

Before releasing software, it's important to ensure each "code" file has a license header, and that
header has an appropriate copyright date. Provided here are two main `bash` scripts to assist
developers in adding and updating the header files.

* `license-prepend-missing.sh` finds "code" files which do not have a license header, and prepends
  the relevant license header to the file.

  *Warning: Matlab uses the header of the file for its automatic help text. For Matlab files,
  license files need to be appended to the file, or removed from the distribution entirely (LIVVkit
  opted for the latter).*

* `license-bump-copyright.sh` will find outdated copyright date statements, and replace them
  (in place) with a current copyright statement.

## Usage

### `license-prepend-missing.sh`

The `license-prepend-missing.sh` script will prepend a license header from a template to files with
these extensions:

| Extension | Template file |
| --------- | ------------- |
| `.py`     | `header-py`  |
| `.sh`     | `header-py`  |
| `.bash`   | `header-py`  |
| `.html`   | `header-html`|
| `.css`    | `header-css` |
| `.js`     | `header-css` |
| `.ncl`    | `header-ncl` |

*Note: Because comment syntax differs between these file types, files mush have one of the above
extensions --- extension-less files, and files without the above extensions are ignored.*

This script work by first sourcing `license-setup.sh`, which sets the important search variables:

* `SOURCE_DIR` sets the top level of the source repository, relative to the current working
  directory. Currently, this is set as `../..` assuming the scripts will be run from this directory
  (`docs/headers`).

* `CURRENT` is a string that should be unique to the license header (e.g., `Copyright (c)`) and is
  used to determine if the header is present in a file.

* `ALWAYS_IGNORE` contains a set of regular expressions to match files which should always be
  ignored by the find command (e.g., `*.git/*`).

* `FILE_IGNORE` contains a set of regular expressions to match the types of files that don't need a
  license header (e.g., `*.png`, `*.md`).

* `PYTHON_IGNORE` contains a set of regular expressions to match python files which don't need a
  license header, or are repackaged and licensed separately.

* `CSS_IGNORE` contains a set of regular expressions to match web (CSS, Javascript, HTML) files
  which don't need a license header, or are repackaged and licensed separately.

The `license-prepend-missing.sh` script then goes through a number of find-and-modify blocks for
the above extensions, and appends a license header to the files missing a header.

**Before** running the `license-prepend-missing.sh` script, it's important to ensure that the script
is finding the correct files.

**First**, run the `license-find-missing.sh` script to verify the found files, which will output a list
of files without a license header.

**Read through this output carefully!**  This find command is generous and not restricted to the
above extensions. Files found with extensions differing from the list above will either need to be
added to the appropriate `*_IGNORE` variables in `license-setup.sh`, or added to the appropriate
find-and-modify block in `license-prepend-missing.sh` (or have a header added to it manually).  *If
this file type's comment syntax differs from those found in the available templates, a new template
will need to be created and a new find-and-modify block added to `license-prepend-missing.sh`.*

Once you're confident the `license-find-missing.sh` script is finding all the appropriate files, and
`license-prepend-missing.sh` has all the required find-and-modify blocks, simply run the
`liscense-prepend-missing.sh` script. You should now see a license header at the top of all your
"code" files.

### `license-bump-copyright.sh`

This script also uses the `*_IGNORE` variables from `license-setup.sh`, but instead of using the
`CURRENT` string to determine if a license header is present, it looks for the `OLD` string, which
is more restrictive, and replaces it (in place) with the `NEW` string.

For example, `CURRENT` is typically set as `Copyright (c)`, which will match files with any copyright
statement, while `OLD` may be set as `Copyright (c) 2015, UT` to match files with the 2015
copyright statement and `NEW` would set to `Copyright (c) 2015,2016, UT` to add 2016 to the
copyright statement. (`OLD` and `NEW` are found in `license-bump-copyright.sh`.)

*Note: Like `license-find-missing.sh`, this script is not restricted to the above file types that
`license-prepend-missing.sh` is.*
