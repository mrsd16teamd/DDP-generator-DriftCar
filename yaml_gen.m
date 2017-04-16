init_control_seq = sprintf('%f,' , U);
init_control_seq = init_control_seq(1:end-1);

outfile = fopen('config.yaml','w');
fprintf(outfile,'init_control_seq: [%s]\n', init_control_seq);
