===csda-dataimputation project===
Copyright 2014 C. P. de Campos (cassiopc@gmail.com). All rights reserved.
This work is licensed under a Creative Commons
Attribution-Noncommercial-Share Alike 3.0 United States License
http://creativecommons.org/licenses/by-nc-sa/3.0/us/

This is the software package to perform data imputation using Bayesian
networks. A Bayesian network model is learned from incomplete data and
then applied to estimate the missing values. The software requires
Matlab, including the BNT package of https://github.com/bayesnet/bnt/.

For technical details, we refer to the paper entitled "Bayesian
network data imputation with application to survival tree analysis",
to appear in the journal Computational Statistics and Data Analysis.
For operational information, please read this file.

It is assumed throughout that you have matlab installed and the BNT
package has been loaded into the matlab path. To run data imputation
using this software, you need to learn a Bayesian network model, which
is accomplished by calling the function 'structureEM' (to learn a
Bayesian network from an incomplete data set using the
expectation-maximization idea). This function also returns the imputed
data set for convenience. A very simple example follows. Further
details on those functions can be found in their corresponding .m
source code files.

% Equivalent sample size for learning the BN
bdeu = 2;
% Load an example file
a=csvread('modelP1r_1500samples.csv');
% each sample must be in a column, so transpose a
a=a';
[n,m] = size(a);
dat=cell(n,m);
% matrix to a cell matrix by hand and create some missing values
% just for testing...
% Furthermore, the data elements must always start in 1,2...
for i = 1:n
  for j = 1:m
    if rand < 0.9
      dat{i,j} = a(i,j)+1;
    else
      dat{i,j} = [];
    end
  end
end

[bnet,ll,datimp,topsort] = structureEM(dat,bdeu,20,20,'EM',1);
% Show the DAG just to visually inspect it
bnet.dag'
% Show the data already imputed
datimp
% To compare to the original data, one needs to reorder the variables
% of the original data first
sum(datimp - a(topsort,:))
