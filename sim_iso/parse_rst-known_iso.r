medium_fim_diag = c(2.88694,8.21349,2.16105,1.25691,2.57742,1.75703,2.10831,2.02828,3.55743);
short_fim_diag = c(2.48775,7.17453,1.18123,0.941488,1.69169,1.20202,1.11135,1.02469,3.10115);

K = length(medium_fim_diag) + 1;

data_dir = "~/bioinfo/data/isoform_mle/sim_iso/";
rst_name = "19869-TCF7-known_iso-mle-med_sho-fixed";
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
isoforms = "Known isoforms";

x11();

plot(0, xlim=c(0, total_cost) * unit_cost, ylim = c(min(vars_fim, vars_actual), max(vars_fim, vars_actual)), type = "n", xlab = "Cost for medium reads", ylab = expression('Average variance of '* (hat(theta[k]))));
title(expression('Average variance of '* (hat(theta[k]))));
text(0.08, 0.12, paste(sep="", "Gene: ", gene, "; Total cost: $", total_cost * unit_cost, "\n", isoforms, "\nMedium and short reads"));
points(unit_cost * c(0:total_cost), vars_fim, col = "red", pch = 1);
lines(unit_cost * c(0:total_cost), vars_fim, col = "red", lty = 1);
points(unit_cost * c(0:total_cost), vars_actual, col = "blue", pch = 2);
lines(unit_cost * c(0:total_cost), vars_actual, col = "blue", lty = 2);
legend(0.03, 0.17, legend = c("Estimation based on FIM", "Simulation result"), col = c("red", "blue"), lty = c(1, 2), pch = c(1, 2));

dev.copy(pdf, paste(sep = "", data_dir, "img/var-known_iso-", gene, ".pdf"));
dev.off();
