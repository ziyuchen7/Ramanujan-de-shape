function [ A ] = Create_Dictionary( Nmax, rowSize, method )
%Create_Dictionary(Nmax,rowSize,method) retruns a dictionary matrix with 
%maximum expected Period = Nmax, and rowSize number of rows.
%method = 'random' or 'Ramanujan' or 'NaturalBasis' or 'DFT'.
%'random' gives a random NPD.
%'Ramanujan' gives the Ramanujan NPD.
%'NaturalBasis' gives the Natural Basis NPD.
%'Farey' gives the Farey NPD.

% The relevant paper is:
% [] S.V. Tenneti and P. P. Vaidyanathan, "Nested Periodic Matrices and Dictionaries:
% New Signal Representations for Period Estimation", IEEE Transactions on Signal 
% Processing, vol.63, no.14, pp.3736-50, July, 2015.

% Copyright (c) 2016, California Institute of Technology. [Based on
% research sponsored by ONR.] All rights reserved.
% Based on work from the lab of P P Vaidyanathan.
% Code developed by Srikanth V. Tenneti.

% No-endorsement clause: Neither the name of the California Institute 
% of Technology (Caltech) nor the names of its contributors may be used to 
% endorse or promote products derived from this software without specific 
% prior written permission.

% Disclaimer: THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND 
% CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
% BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



A=[];

for N=1:Nmax    
    
    if(strcmp(method,'Farey')==0)
        
    if(strcmp(method,'Ramanujan')==1)
        c1 = zeros(N,1);
        k_orig = 1:N; k = k_orig(gcd(k_orig,N)==1);

        for n = 0:(N-1)
         for a = k
          c1(n+1) = c1(n+1) + exp(1j*2*pi*a*(n)/N);
         end
        end

        c1=real(c1);
    elseif(strcmp(method,'NaturalBasis')==1)
        c1 = zeros(N,1);c1(1) = 1;
    elseif(strcmp(method,'random')==1)
        c1 = randn(N,1);
    end
        
    k_orig=1:N;k=k_orig(gcd(k_orig,N)==1);CN_colSize = size(k,2);
    CN=[];
    
    for j = 1:CN_colSize
        CN = cat(2,CN,circshift(c1,(j-1)));
    end
    
     else
        A_dft = dftmtx(N);
        a = 0:N-1;a(1)=N; a=N./gcd(a,N); 
        I = 1:N;
        I = I(a(I)==N);
        CN=A_dft(:,I);
    end
    
    CNA = repmat(CN,floor(rowSize/N),1);
    CN_cutoff = CN(1:rem(rowSize,N),:);
    CNA =cat(1,CNA,CN_cutoff);
    
   
   
    A=cat(2,A,CNA);
end

A = round(A);


end

