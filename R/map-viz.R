#' Coerce a raw core-metadata data frame into a `sample_registry`
#'
#' Maps arbitrary column names from a raw input file (e.g. the
#' `core_metadata.tsv` shape used across sedaDNA projects, one row per
#' subsample) onto the standardized [new_sample_registry()] schema. Column
#' name arguments default to the names used in that convention; override
#' them for differently-named source files.
#'
#' Only `core_id`, `depth_cm`, `latitude`, and `longitude` are required --
#' `core_id` defaults to reading from a column named `"sample_source"`
#' (the `core_metadata.tsv` convention's true physical-core grouping key;
#' `core_ID`/`sample_ID` are per-subsample composite identifiers like
#' `"GeoB9064-1-251-A"`, not this). `water_depth_m` and `age_ka` are
#' filled `NA` if their column is absent; `sample_id` is auto-generated
#' (`"<core_id>_<row number>"`) if absent; `core_label` is left out of
#' the result entirely if absent (rather than filled `NA`) since
#' [collapse_to_core_registry()] falls back to `core_id` itself when
#' there's no separate label.
#'
#' @param x A raw data frame, one row per subsample.
#' @param core_id,depth_cm,latitude,longitude Names of the columns in `x`
#'   holding each required field.
#' @param water_depth_m,sample_id,age_ka Names of the columns in `x`
#'   holding each optional field; used if present, otherwise filled `NA`
#'   (or auto-generated, for `sample_id`). Set to `NULL` to always skip
#'   the lookup.
#' @param core_label Name of an optional column holding a human-readable
#'   core label (e.g. `"GeoB1706-2; 2003"`, often including the retrieval
#'   year); kept as a passthrough column and used by
#'   [collapse_to_core_registry()]. Set to `NULL` if absent.
#' @return A `sample_registry` tibble.
#' @export
as_sample_registry <- function(x,
                                core_id = "sample_source",
                                depth_cm = "depth",
                                latitude = "latitude",
                                longitude = "longitude",
                                water_depth_m = "water_depth",
                                sample_id = "sample_ID",
                                age_ka = "age",
                                core_label = "core_labels") {
  required_cols <- c(
    core_id = core_id, depth_cm = depth_cm, latitude = latitude, longitude = longitude
  )
  missing_cols <- setdiff(unname(required_cols), names(x))
  if (length(missing_cols) > 0) {
    rlang::abort(sprintf(
      "Column(s) not found in `x`: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  x <- tibble::as_tibble(x)

  core_id_vals <- as.character(x[[required_cols[["core_id"]]]])

  out <- tibble::tibble(
    core_id = core_id_vals,
    depth_cm = as.numeric(x[[required_cols[["depth_cm"]]]]),
    latitude = as.numeric(x[[required_cols[["latitude"]]]]),
    longitude = as.numeric(x[[required_cols[["longitude"]]]]),
    water_depth_m = if (!is.null(water_depth_m) && water_depth_m %in% names(x)) {
      as.numeric(x[[water_depth_m]])
    } else {
      NA_real_
    },
    sample_id = if (!is.null(sample_id) && sample_id %in% names(x)) {
      as.character(x[[sample_id]])
    } else {
      paste(core_id_vals, seq_len(nrow(x)), sep = "_")
    },
    age_ka = if (!is.null(age_ka) && age_ka %in% names(x)) {
      as.numeric(x[[age_ka]])
    } else {
      NA_real_
    }
  )

  if (!is.null(core_label) && core_label %in% names(x)) {
    out$core_label <- x[[core_label]]
  }

  new_sample_registry(out)
}

#' Collapse a `sample_registry` to one row per core
#'
#' Generalizes the `separate(core_labels, ...) |> distinct()` pattern used
#' to derive a core-level table for mapping and retrieval: assumes
#' latitude/longitude/water depth are constant within a core (true for all
#' piston/gravity core sampling in this domain) and takes the first
#' non-missing value per core.
#'
#' @param x A `sample_registry`, as produced by [as_sample_registry()].
#' @return A `core_registry` tibble, one row per `core_id`.
#' @export
collapse_to_core_registry <- function(x) {
  validate_sample_registry(x)

  core_label <- if ("core_label" %in% names(x)) x$core_label else x$core_id
  x <- dplyr::mutate(x, .core_label = core_label)

  out <- x |>
    dplyr::group_by(.data$core_id) |>
    dplyr::summarise(
      core_label = dplyr::first(.data$.core_label),
      latitude = dplyr::first(.data$latitude),
      longitude = dplyr::first(.data$longitude),
      water_depth_m = dplyr::first(.data$water_depth_m),
      n_samples = dplyr::n(),
      depth_min_cm = min(.data$depth_cm, na.rm = TRUE),
      depth_max_cm = max(.data$depth_cm, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      source_program = classify_core_source(.data$core_id),
      pangaea_doi = NA_character_,
      n_samples = as.integer(.data$n_samples)
    )

  new_core_registry(out)
}

#' Plot cores on a static world map
#'
#' Ports the `rnaturalearth` + `ggrepel` map pattern used across sedaDNA
#' projects into a single reusable function.
#'
#' @param x A `core_registry`, as produced by [collapse_to_core_registry()].
#' @param xlim,ylim Optional longitude/latitude limits; default to the data
#'   extent padded by 10 degrees.
#' @param label Logical; add repelled core-ID labels?
#' @param point_size Point size for core markers.
#' @return A `ggplot` object.
#' @export
plot_core_map <- function(x, xlim = NULL, ylim = NULL, label = TRUE, point_size = 5) {
  validate_core_registry(x)

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  if (is.null(xlim)) xlim <- range(x$longitude, na.rm = TRUE) + c(-10, 10)
  if (is.null(ylim)) ylim <- range(x$latitude, na.rm = TRUE) + c(-10, 10)

  p <- ggplot2::ggplot(world) +
    ggplot2::geom_sf(fill = "#272733") +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::geom_point(
      data = x,
      ggplot2::aes(x = .data$longitude, y = .data$latitude, fill = .data$core_id),
      size = point_size, shape = 21, alpha = 0.85
    ) +
    ggplot2::scale_fill_viridis_d(name = "Core") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = NULL)

  if (label) {
    p <- p + ggrepel::geom_label_repel(
      data = x,
      ggplot2::aes(x = .data$longitude, y = .data$latitude, label = .data$core_id),
      color = "black", fill = "white", size = 2.5, box.padding = 0.4
    )
  }

  p
}

#' Marker colors for [leaflet_core_map()] / [plot_core_map()]
#'
#' Both chosen to stand out against the blue/teal ocean basemap tiles,
#' which a default ggplot2/leaflet blue does not.
#' @noRd
map_colors <- c(input = "#e6550d", neighbour = "#984ea3")

#' Plot cores on an interactive leaflet map
#'
#' @param x A `core_registry`, as produced by [collapse_to_core_registry()]
#'   -- the input/query cores, plotted in orange.
#' @param neighbours An optional `neighbour_registry`, as produced by
#'   [find_neighbours()], plotted as a second layer in purple so
#'   neighbouring cores are visually distinguishable from the input set.
#'   Rows with missing coordinates are dropped.
#' @return A `leaflet` htmlwidget.
#' @export
leaflet_core_map <- function(x, neighbours = NULL) {
  validate_core_registry(x)

  map <- leaflet::leaflet() |>
    leaflet::addProviderTiles(leaflet::providers$Esri.OceanBasemap) |>
    leaflet::addCircleMarkers(
      data = x,
      lng = ~longitude, lat = ~latitude,
      radius = 6, stroke = TRUE, color = map_colors[["input"]], fillOpacity = 0.8,
      label = ~core_id,
      popup = ~sprintf(
        "<b>%s</b><br/>Water depth: %s m<br/>Samples: %s<br/>Program: %s",
        core_id, water_depth_m, n_samples, source_program
      ),
      layerId = ~core_id,
      group = "Input cores"
    )

  if (!is.null(neighbours)) {
    validate_neighbour_registry(neighbours)
    neighbours <- neighbours[!is.na(neighbours$latitude) & !is.na(neighbours$longitude), ]
  }

  if (!is.null(neighbours) && nrow(neighbours) > 0) {
    map <- map |>
      leaflet::addCircleMarkers(
        data = neighbours,
        lng = ~longitude, lat = ~latitude,
        radius = 5, stroke = TRUE, color = map_colors[["neighbour"]], fillOpacity = 0.7,
        label = ~neighbour_core_id,
        popup = ~sprintf(
          "<b>%s</b><br/>Distance: %.1f km<br/>Program: %s<br/>DOI: %s",
          neighbour_core_id, distance_km, source_program, neighbour_doi
        ),
        layerId = ~paste0("neighbour-", neighbour_doi),
        group = "Neighbouring cores"
      ) |>
      leaflet::addLegend(
        position = "bottomright",
        colors = unname(map_colors),
        labels = c("Input cores", "Neighbouring cores")
      )
  }

  map
}
