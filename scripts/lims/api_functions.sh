#!/bin/bash

# Requires the env vars LIMS_API_URL and LIMS_API_TOKEN to be set
# Requires `jq` to be on your PATH

lims_get_by_url() {
  /usr/bin/curl \
    "$1" \
    -H "Authorization: Token $LIMS_API_TOKEN" \
    -k \
    2>/dev/null
}

lims_get () {
  lims_get_by_url "$LIMS_API_URL/$1"
}

lims_get_all () {
  response=$(lims_get $1)
  jq '.results[]' <<< $response
  next=$(jq -r .next <<< $response)
  while [ "$next" != "null" ] ; do
    response=$(lims_get_by_url $next )
    jq '.results[]' <<< $response
    next=$(jq -r .next <<< $response )
  done
}
