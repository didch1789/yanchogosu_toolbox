function [S, Q] = multslice_community_detection(A, 
% For studying!!
% combine Genlouvain with BCT!
A = fmri_gen_pseudodat_ycgosu(1000, 100, 'noisetype', 'gaussian');
out = struct;
for i = 1:5
    out.(sprintf('A%d', i)) =corr(A(i*200-199:i*200, :));
end
out2 = struct2cell(out);

gamma = 1;
omega = 1;

N=length(out2{1});
T=length(out2);
B = sparse(N*T, N*T); 
twomu=0;

for s=1:T
    k=sum(out2{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=out2{s}-gamma*k'*k/twom;
end
twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B);
Q = Q/twomu
S = reshape(S,N,T);