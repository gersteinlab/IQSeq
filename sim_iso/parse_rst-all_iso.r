medium_fim_diag = c(2.91649,2.89989,2.39377,11.27,2.88871,2.27176,2.90086,2.88313,2.35979,11.5711,2.87028,2.23267,2.64738,2.62227,2.03648,11.7312,2.59602,1.87866,2.14864,2.11428,1.50965,9.49094,2.06822,1.31648,2.14714,2.10739,1.41921,11.0271,2.0533,1.20128,1.59602,1.54483,0.832266,7.73876,1.46619,0.570323,2.69356,2.66841,2.0686,12.2844,2.64275,1.90985,2.18417,2.14947,1.53014,9.94056,2.10334,1.33486,2.18672,2.14652,1.44047,11.6612,2.09229,1.21985,1.59463,1.54208,0.810858,7.96891,1.46127,0.542738,2.53586,2.50625,1.87435,12.1621,2.47218,1.69729,1.99074,1.95058,1.29757,9.40735,1.89386,1.08013,2.20599,2.16558,1.45088,11.982,2.11129,1.22898,1.59395,1.54074,0.800489,8.08329,1.45888,0.529374,3.57035,3.56091,2.99512,15.7512,3.56874,2.89639,1.06826,0.980921,0.775742,42.6546,0.545459);
short_fim_diag = c(2.31116,2.30352,2.33633,10.3273,2.28792,2.17164,2.15041,2.14012,2.17483,10.0474,2.11861,1.9968,1.40955,1.38673,1.42938,8.10998,1.33713,1.19204,1.19162,1.16479,1.21004,7.47211,1.10619,0.953494,1.00357,0.970685,1.02264,7.84456,0.898472,0.726648,0.734185,0.695889,0.75125,6.808,0.611451,0.429655,1.4135,1.39026,1.43381,8.35238,1.33973,1.19208,1.19092,1.16355,1.20979,7.69045,1.10373,0.948231,0.988741,0.954947,1.00818,8.0495,0.880676,0.705246,0.712581,0.673166,0.729935,6.95559,0.586203,0.400434,1.20641,1.17889,1.22576,7.92919,1.11872,0.961456,0.965045,0.932892,0.982717,7.12504,0.862306,0.696397,0.981576,0.947339,1.00119,8.15051,0.87207,0.694898,0.702137,0.66218,0.719628,7.02795,0.573989,0.3863,2.73347,2.73056,2.76485,13.2831,2.72571,2.60579,0.533936,0.465033,0.566917,26.8963,0.310591);

K = length(medium_fim_diag) + 1;

medium_read_length = 250;
medium_bp_cost = 7E-5;
short_read_length = 30;
short_bp_cost = 7E-6;
epsilon = 1E-9;

data_dir = "~/bioinfo/data/isoform_mle/sim_iso/";
rst_name = "19869-TCF7-all_iso-mle-med_sho-fixed";
unit_cost = 0.01;
total_cost = 20;

rst_file = paste(sep = "", data_dir, rst_name, "-unit_", unit_cost, "-total_", total_cost, ".rst");
raw_rst = read.table(rst_file, header = FALSE, col.names = c("medium_read_cost", "short_read_cost", "var"));

vars_actual = c();
vars_fim = c();
for (medium_cost in 0:total_cost) {
	short_cost = total_cost - medium_cost;
	var_actual = subset(raw_rst, medium_read_cost == medium_cost)$var[1];
	vars_actual = append(vars_actual, var_actual);

	num_medium_reads = medium_cost * unit_cost / medium_bp_cost / medium_read_length;
	num_short_reads = short_cost * unit_cost / short_bp_cost / short_read_length;
	comb_diag = num_medium_reads * medium_fim_diag + num_short_reads * short_fim_diag;
	vars_fim = append(vars_fim, sqrt(sum(1 / comb_diag) / (K - 1)));
}

gene = "TCF7";
isoforms = "All possible isoforms";

x11();

plot(0, xlim=c(0, total_cost) * unit_cost, ylim = c(min(vars_fim, vars_actual), max(vars_fim, vars_actual)), type = "n", xlab = "Cost for medium reads", ylab = expression('Average variance of '* (hat(theta[k]))));
title(expression('Average variance of '* (hat(theta[k]))));
text(0.08, 0.12, paste(sep="", "Gene: ", gene, "; Total cost: $", total_cost * unit_cost, "\n", isoforms, "\nMedium and short reads"));
points(unit_cost * c(0:total_cost), vars_fim, col = "red", pch = 1);
lines(unit_cost * c(0:total_cost), vars_fim, col = "red", lty = 1);
points(unit_cost * c(0:total_cost), vars_actual, col = "blue", pch = 2);
lines(unit_cost * c(0:total_cost), vars_actual, col = "blue", lty = 2);
legend(0.03, 0.2, legend = c("Estimation based on FIM", "Simulation result"), col = c("red", "blue"), lty = c(1, 2), pch = c(1, 2));

dev.copy(pdf, paste(sep = "", data_dir, "img/var-all_iso-", gene, ".pdf"));
dev.off();
