#!/bin/bash

# This file is a library intended for use with http://sentry.io
# Example usage in your script:
: <<USAGE_EXAMPLE
#!/bin/bash
source sentry-lib.bash
the_rest_of_your_script
USAGE_EXAMPLE

# Configuration:
# One of two ways:
# 1) Create $HOME/.sentryclirc
# 2) Set ENV variable
# These require a DSN which can be found at https://sentry.io
# https://sentry.io/settings/[your_org]/projects/[your_project]/keys/
: <<RC_EXAMPLE
[auth]
token=redacted
dsn=redacted

[defaults]
org=your-org-goes-here
project=your-project-goes-here
RC_EXAMPLE

: <<VARS_EXAMPLE
export SENTRY_DSN=<your-dsn-goes-here>
VARS_EXAMPLE

# Guard against direct execution
(return 0 2>/dev/null) || (
    echo >&2 "Error: You should source this, like 'source \"$0\"', not run this directly. See comments for usage."
    exit 1
)

if ! command -v sentry-cli >/dev/null ; then
    echo >&2 "WARNING: sentry-cli not available on your PATH"
elif ! [[ -s "$HOME/.sentryclirc" ]] && [[ "$SENTRY_DSN" == "" ]]  ; then
    echo >&2 "WARNING: sentry-cli improperly configured"
else
    eval "$(sentry-cli bash-hook)"
fi
