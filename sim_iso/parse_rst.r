rst_names = c(
"19869-TCF7-200901302216-med_sho",
"19869-TCF7-200901311021-med_sho",
"19869-TCF7-200901312218-med_sho-fixed",
"19869-TCF7-200901312207-med_sho-fixed",
"19869-TCF7-200902012002-med_sho-fixed2",
"19869-TCF7-200902012005-med_sho-fixed2");

rst_descs = c(
"TCF7: random probabilities",
"TCF7: random probabilities",
"TCF7: fixed probabilities: equal",
"TCF7: fixed probabilities: equal",
"TCF7: fixed probabilities: uneven: 0.8, 0.1, ...",
"TCF7: fixed probabilities: uneven: 0.8, 0.1, ...");

img_file_names = c(
"TCF7-rand_prob",
"TCF7-fixed_prob_eq",
"TCF7-fixed_prob_unev"
);

unit_costs = c(
0.01,
0.02,
0.01,
0.02,
0.01,
0.02);

total_costs = c(
20,
20,
20,
20,
20,
20);

err_idx = 3;
err_type = "known MAE";

#err_idx = 6;
#err_type = "all RMSE";

data_dir = "~/bioinfo/data/isoform_mle/sim_iso/";

rst_idx = 1:6;

for (i in rst_idx) {
	if (i %% 2 == 1) {
		x11();
		par(mfrow=c(2,1));
	}
	rst_name = rst_names[i];
	unit_cost = unit_costs[i];
	total_cost = total_costs[i];
	rst_desc = rst_descs[i];

	rst_file = paste(sep = "", data_dir, rst_name, "-unit_", unit_cost, "-total_", total_cost, ".rst");

	raw_rst = read.table(rst_file, header = FALSE, col.names = c("medium_read_cost", "short_read_cost", "known_mae", "known_rmse", "all_mae", "all_rmse"));

	epsilon = 1E-9;
	cost_counts = 0:total_cost;
	costs = cost_counts * unit_cost;
	num_costs = length(costs);

	errs = list();
	mean_errs = list();
	labels = c();
	for (cost_med in cost_counts) {
		cost_sho = total_cost - cost_med;
		errs[[cost_med + 1]] = raw_rst[which(raw_rst[, 1] == cost_med), err_idx];
		mean_errs[[cost_med + 1]] = mean(errs[[cost_med + 1]]);
		labels = c(labels, paste(sep = "", "m:$", cost_med*unit_cost,
			"\ns:$", cost_sho*unit_cost));
	}

	best_idx = which(unlist(mean_errs) == min(unlist(mean_errs)))[1];

	boxplot(x = errs, names = labels, cex.axis = 0.7);
	boxplot(x = mean_errs, names = NA, add = TRUE, border = "red", yaxt = "n");
	boxplot(x = mean_errs[[best_idx]], at = best_idx, names = NA, add = TRUE, border = "blue", yaxt = "n");
	legend(1, max(unlist(errs)) * 0.95, col = c("red", "blue"), lty = c(1, 1), legend = c("Mean estimation error", paste("lowest mean estimation error:", mean_errs[[best_idx]])));
	title(main = paste(sep = "", "Average estimation error: ", err_type , " (w/ total cost $", total_cost * unit_cost, ")\n", rst_desc), sub = "$0.2 medium reads: ~1x, $0.2 short reads: ~10x");
	# par(ask=TRUE);
	if (i %% 2 == 0) {
		dev.copy(pdf, paste(sep = "", data_dir, "img/", img_file_names[i %/% 2], ".pdf"));
		dev.off();
	}
}
