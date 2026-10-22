# Independently authored, fixed-data executable plumbing check. No benchmark.
# Usage: Rscript scripts/run_stage1_smoke.R PACKAGE_ROOT OUTPUT_DIR
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
package_root <- normalizePath(args[[1]])
output <- args[[2]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
devtools::load_all(package_root, quiet = TRUE)
internal <- function(name) getFromNamespace(name, "rajiveplus")
hash <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)
gram <- function(q) max(abs(crossprod(q) - diag(ncol(q))), 0)
basis <- internal(".ablation_basis")
write_json <- function(x, name) jsonlite::write_json(
  x, file.path(output, name), pretty = TRUE, auto_unbox = TRUE,
  digits = 16, null = "null", na = "null"
)
seeds <- list(data = 1729L, block_fit = 8101L, aggregation = 8102L, individual = 8103L)
manifest_path <- file.path("config", "stage1_ablation.json")
manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
spec <- manifest$seed_manifest$smoke
stopifnot(identical(as.integer(spec$seeds), seeds$data),
          spec$stage_seeds$block_fit == seeds$block_fit,
          spec$stage_seeds$aggregation == seeds$aggregation,
          spec$stage_seeds$individual_selection == seeds$individual,
          spec$dimensions$n == 36, all(spec$dimensions$p == c(12, 10, 8)),
          all(spec$latent_generator_arguments$joint_strength == c(4, 3)),
          spec$latent_generator_arguments$individual_strength == 2,
          spec$latent_generator_arguments$snr == 4)
rank_tol <- function(x) max(dim(x)) * .Machine$double.eps *
  max(svd(x, nu = 0L, nv = 0L)$d, 0)
sim <- internal(".simulate_identifiable_data")(
  K = 3L, n = 36L, pks = c(12L, 10L, 8L), rank_joint = 2L,
  rank_individual = c(1L, 1L, 1L), joint_strength = c(4, 3),
  individual_strength = 2, snr = 4, seed = seeds$data
)
certificate <- internal(".certify_simulated_truth")(sim)
stopifnot(certificate$status == "PASS")
write_json(certificate, "smoke_truth_certificate.json")
cache <- internal(".ablation_prepare")(sim$noisy, c(3L, 3L, 3L), seeds$block_fit)
data_hash <- hash(sim$noisy)
fit_hash <- hash(cache$block_fits)
rows <- list()
details <- list()
baseline_result <- NULL
for (R in c("raw", "orthonormal")) for (S in c("robust", "classical")) {
  geometry <- tryCatch(internal(".ablation_geometry")(
    cache, 2L, R, S, seeds$aggregation
  ), error = identity)
  selected <- if (!inherits(geometry, "error")) tryCatch(
    internal(".ablation_select_individual")(cache, geometry, seeds$individual),
    error = identity
  ) else geometry
  for (I in c("original", "constrained")) {
    id <- paste(R, S, I, sep = "__")
    failure <- if (inherits(selected, "error")) selected else NULL
    result <- if (is.null(failure)) tryCatch(internal(".ablation_finish")(
      cache, geometry, I, seeds$individual, selected
    ), error = identity) else failure
    if (inherits(result, "error")) {
      rows[[id]] <- data.frame(
        arm = id, R = R, S = S, I = I, status = "FAIL",
        pre_filter_rank = NA_integer_, post_filter_rank = NA_integer_,
        joint_gram_error = NA_real_, block_gram_error = NA_real_,
        maximum_span_leakage = NA_real_, maximum_rank_drop = NA_integer_,
        maximum_closure = NA_real_
      )
      details[[id]] <- list(error = conditionMessage(result))
      next
    }
    q <- geometry$joint_scores
    rows[[id]] <- data.frame(
      arm = id, R = R, S = S, I = I, status = "SMOKE_PASS_NOT_CALIBRATED",
      pre_filter_rank = geometry$pre_filter_rank,
      post_filter_rank = geometry$post_filter_rank,
      joint_gram_error = gram(q),
      block_gram_error = max(vapply(geometry$block_scores, gram, numeric(1L))),
      maximum_span_leakage = max(vapply(result$blocks, `[[`, numeric(1L),
                                        "individual_fitted_span_leakage")),
      maximum_rank_drop = max(vapply(result$blocks, `[[`, integer(1L),
                                     "individual_rank_drop")),
      maximum_closure = max(vapply(result$blocks, `[[`, numeric(1L),
                                   "reconstruction_relative"))
    )
    details[[id]] <- list(
      data_sha256 = data_hash, block_fit_sha256 = fit_hash,
      selected_individual_sha256 = hash(selected),
      stack_sha256 = hash(geometry$stacked_scores),
      dropped_directions = geometry$identif_dropped,
      retained_directions = setdiff(seq_len(geometry$pre_filter_rank),
                                   geometry$identif_dropped),
      joint_gram_pre_filter = gram(geometry$pre_filter_scores),
      joint_gram_post_filter = gram(geometry$joint_scores),
      representation_algorithm = if (R == "raw") "unchanged robust coordinates"
        else "base::svd left singular vectors; native LAPACK order and signs",
      representation_rank_tolerance_policy = "max(dim)*machine_epsilon*largest_singular_value",
      block_representation_tolerances = vapply(seq_along(cache$blocks), function(k) {
        rank_tol(cache$block_fits[[k]]$u[, seq_len(cache$initial_ranks[k]), drop = FALSE])
      }, numeric(1L)),
      final_representation_tolerance = rank_tol(geometry$raw_selected_scores),
      representation_rank_drop = geometry$final_representation_rank_drop,
      requested_joint_rank = 2L, achieved_joint_rank = ncol(q),
      individual_ranks_before = vapply(selected, `[[`, integer(1L), "rank"),
      individual_ranks_after = vapply(result$blocks, function(b) b$individual$rank,
                                      integer(1L)),
      projector_interpretation = geometry$projector_interpretation
    )
    if (R == "raw" && S == "robust" && I == "original") baseline_result <- result
  }
}
table <- do.call(rbind, rows)
rownames(table) <- NULL
write.csv(table, file.path(output, "eight_arm_smoke.csv"), row.names = FALSE)
write_json(details, "smoke_pairing_manifest.json")
# Public baseline invocation has identical preprocessing (none), robust block
# construction, fixed initial/target ranks and l2 identifiability.
public <- Rajive(sim$noisy, c(3L, 3L, 3L), joint_rank = 2L, seed = seeds$block_fit,
                 identifiability_norm = "l2", num_cores = 1L)
stopifnot(!is.null(baseline_result))
baseline_checks <- list(
  joint_scores = isTRUE(all.equal(public$joint_scores,
    baseline_result$geometry$joint_scores, tolerance = 1e-12)),
  post_filter_rank = identical(public$joint_rank,
    baseline_result$geometry$post_filter_rank)
)
for (k in seq_along(sim$noisy)) for (role in c("joint", "individual", "residual")) {
  actual <- internal(".get_component_matrix")(public, k, role)
  expected <- baseline_result$blocks[[k]][[role]]
  if (role != "residual") expected <- expected$full
  baseline_checks[[paste(k, role, sep = "_")]] <-
    isTRUE(all.equal(unname(actual), unname(expected), tolerance = 1e-12))
  baseline_checks[[paste(k, role, "public_dimnames", sep = "_")]] <-
    identical(dimnames(actual), dimnames(sim$noisy[[k]]))
}
passed <- all(table$status == "SMOKE_PASS_NOT_CALIBRATED") &&
  all(unlist(baseline_checks)) &&
  all(table$joint_gram_error[table$R == "orthonormal"] <= 1e-10) &&
  all(table$block_gram_error[table$R == "orthonormal"] <= 1e-10) &&
  all(table$maximum_span_leakage[table$I == "constrained"] <= 1e-8) &&
  all(table$maximum_closure <= 1e-12)
write_json(list(
  status = if (passed) "PASS_CORRECTNESS_ONLY" else "FAIL",
  scope = "one synthetic fixed-data eight-arm smoke; no performance estimate",
  seeds = seeds, n = 36L, p = c(12L, 10L, 8L),
  target_joint_rank = 2L, initial_ranks = c(3L, 3L, 3L),
  generator = sim$metadata, preprocessing = cache$preprocessing,
  truth_certificate_status = certificate$status,
  input_sha256 = data_hash, block_fit_sha256 = fit_hash,
  protocol_config_sha256 = digest::digest(file = manifest_path, algo = "sha256"),
  package_commit = system2("git", c("-C", package_root, "rev-parse", "HEAD"), stdout = TRUE),
  hash_encoding = "SHA256 of R serialization via digest; matched R environment required",
  baseline_equivalence = baseline_checks,
  arms = nrow(table), software_readiness = "NOT_READY",
  automatic_rank_layer = "BLOCKED_PENDING_MATCHED_NULL",
  session_info = capture.output(sessionInfo())
), "smoke_status.json")
print(table)
stopifnot(passed)
