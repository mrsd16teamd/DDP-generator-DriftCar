T_horizon = 150;

xDes = [8,0.3,0,0,0,0];
X_des = sprintf('%f,' , xDes);
X_des = X_des(1:end-1);

infile = 'U_standstill.txt'
U = load(infile);
U = U';

init_control_seq = sprintf('%f,' , U);
init_control_seq = init_control_seq(1:end-1);

outfile = fopen('config.yaml','w');
fprintf(outfile,'T_horizon: %d\n', T_horizon);
fprintf(outfile,'X_des: [%s]\n', X_des);
fprintf(outfile,'init_control_seq: [%s]\n', init_control_seq);
