library("vcr") # *Required* as vcr is set up on loading
invisible(vcr::vcr_configure(
  dir = vcr::vcr_test_path("fixtures"),
  # yaml (the default) is unreliable here: PANGAEA responses embed long
  # JSON-LD strings with \uXXXX escapes, and vcr/yaml's line-wrapping of
  # long double-quoted scalars can split a \uXXXX sequence across lines,
  # producing a cassette file that fails to parse back in a fresh session
  # even though it round-trips fine within the recording session.
  serialize_with = "json"
))
