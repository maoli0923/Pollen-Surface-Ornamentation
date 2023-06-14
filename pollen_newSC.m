%compute k-dim features based on subgraph centrality
%input: h is binary image (black and white). 
%output: g is k-dim feature

function [g]=pollen_newSC(hh)
G=imageGraph(size(hh),4);
E=G.Edges;
f=hh(E.EndNodes(:,1))+hh(E.EndNodes(:,2));
a=find(f==1);
wt=ones(size(f));
wt(a)=0.01;
G.Edges.Weight=wt;
nn = numnodes(G);
[s,t] = findedge(G);
W = sparse(s,t,G.Edges.Weight,nn,nn);
W=W+W';
%subgraph centrality
[phi lambda]=eig(full(W));
lambda = diag(lambda);
[Val Ind]=sort(sum(phi.^2*exp(lambda),2),'descend');

x=0:5:100;% this is 20-dim. e.g. use 0:2:100 for 50-dim
g=x;
g(1)=0;
for i=2:length(g)
    t=round(x(i)*size(W,1)/100);
    sub=Ind(1:t);
    H = subgraph(G,sub);
    bins = conncomp(H);
    g(i)=max(bins);
end
g=g(2:end);
end