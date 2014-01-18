#[LOG 2009-06-09 19:13:49 INFO] Known isoform probabilities: 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
#[LOG 2009-06-09 19:17:16 INFO] BruteforceFIM(Read, PossibleIsoforms)
#[LOG 2009-06-09 19:17:16 INFO] # ofim calls: 26902
#[LOG 2009-06-09 19:19:09 INFO] FastFIM(Read, PossibleIsoforms)
#[LOG 2009-06-09 19:19:09 INFO] # ofim calls: 169
#[LOG 2009-06-09 19:19:09 INFO] BruteforceFIM = FastFIM
#[LOG 2009-06-09 19:19:35 INFO] FasterFIM(Read, PossibleIsoforms)
#[LOG 2009-06-09 19:19:35 INFO] # ofim calls: 46

gene = "TCF7";
read = "30bp short read";
isoforms = "All possible isoforms";
data_dir = "~/bioinfo/data/isoform_mle/sim_iso/"

num_ofim_calls = c(26902, 169, 46);
exec_times = c((3*60+27), (1*60+53) / 100, 26 / 53);
algs = c("Brute-force FIM", "FastFIM", "FasterFIM");

ofim_speedups = num_ofim_calls[1] / num_ofim_calls;
exec_speedups = exec_times[1] / exec_times;

x11();
boxplot(as.list(ofim_speedups), boxwex = 0.2, names = algs, border = "red", log = "y", sub = paste(sep = "", "Gene: ", gene, "; Read: ", read, "\n", isoforms));
title("Speedups in computing FIM");
lines(1:3, ofim_speedups, col = "red");
lines(1:3, exec_speedups, col = "blue", lty = 2);
boxplot(as.list(exec_speedups), boxwex = 0.2, names = rep("", 3), border = "blue", add = TRUE, log = "y");
legend(1.5, 2.5, legend = c("Theoretical speedup", "Actual speedup in a typical run"), col = c("red", "blue"), lty = c(1, 2));
dev.copy(pdf, paste(sep = "", data_dir, "img/fim-speedup-", gene, ".pdf"));
dev.off();
