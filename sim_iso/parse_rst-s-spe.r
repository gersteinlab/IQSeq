data_dir = "~/bioinfo/data/isoform_mle/sim_iso/";
isoforms = "Known isoforms";

unit_cost = 0.01;
total_cost = 20;

#rst_name = "19869-TCF7-known_iso-mle-20090623-2331";
#rst_name = "19869-TCF7-known_iso-mle-20090624-0027";
#rst_name = "19869-TCF7-known_iso-mle-20090624-0054";
#gene = "TCF7";

#rst_name = "11314-BPTF-known_iso-mle-20090623-2330";
#rst_name = "11314-BPTF-known_iso-mle-20090624-0027";
#rst_name = "11314-BPTF-known_iso-mle-20090624-0054";
rst_name = "11314-BPTF-known_iso-mle-20090624-0207"; unit_cost = 0.015;
#rst_name = "11314-BPTF-known_iso-mle-20090624-0225"; unit_cost = 0.03;
gene = "BPTF";

#rst_name = "10241-WSCD1-known_iso-mle-20090623-2329";
#rst_name = "10241-WSCD1-known_iso-mle-20090624-0054";
#gene = "WSCD1";

rst_file = paste(sep = "", data_dir, rst_name, ".rst");
raw_rst = read.table(rst_file, header = FALSE, col.names = c("medium_read_cost", "short_read_cost", "spe_read_cost", "var"));

vars_actual = c();
for (spe_cost in 0:total_cost) {
	var_actual = subset(raw_rst, spe_read_cost == spe_cost)$var[1];
	vars_actual = append(vars_actual, var_actual);
}

x11();

plot(0, xlim=c(0, total_cost) * unit_cost, ylim = c(min(vars_actual), max(vars_actual)), type = "n", xlab = "Cost for short PE reads \(\$\)", ylab = expression('Average variance of '* (hat(theta[k]))));
title(expression('Average variance of '* (hat(theta[k]))));
text(unit_cost * 8, 0.8*(max(vars_actual)-min(vars_actual))+min(vars_actual), paste(sep="", "Gene: ", gene, "; Total cost: $", total_cost * unit_cost, "\n", isoforms, "\nshort and short PE reads"));
points(unit_cost * c(0:total_cost), vars_actual, col = "blue", pch = 2);
lines(unit_cost * c(0:total_cost), vars_actual, col = "blue", lty = 2);

dev.copy(pdf, paste(sep = "", data_dir, "img/var-known_iso-s-spe-", gene, "-unit_cost-", unit_cost, ".pdf"));
dev.off();
