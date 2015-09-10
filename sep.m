% SEP goes once through the scrambled mixed signals, x
% (which is of length P), in batch blocks of size B, adjusting weights,
% w, at the end of each block.
sweep=sweep+1; t=1;
noblocks=fix(P/B);
BI=B*Id;
for t=t:B:t-1+noblocks*B,
u=w*x(:,t:t+B-1);
w=w+L*(BI+(1-2*(1./(1+exp(-u))))*u')*w;
end;
sepout
