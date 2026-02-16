function [J,BF,e,e0]=seq_glitch_det(z,fran,fdran, J, BF, m)
%
% N=length(z);
consec = 0;
if(nargin < 4)
    J=[];
    m=0;
    BF=[];
end
N=length(z);
p=2; %safeguard at the edges TODO safeguarding
while true
	m=m+1;
	% zero hypothesis
	[~,~,alpha0,~,evidence0,lik] = HMM_Pulse_3dcol(z,fran,fdran,J);
    E0 = evidence0(end)
	e0(m)=E0;
	E=E0*ones(1,N);
	%look for biggest glitch
    tic
	for k=p:N-p+1
        evidence = [];
        if(consec)
            [~,evidence] = HMM_Pulse_3dcol_seq(z,fran,fdran,[J k k+1],alpha0,evidence0,lik,k);
        else
            [~,evidence] = HMM_Pulse_3dcol_seq(z,fran,fdran,[J k],alpha0,evidence0,lik,k);
        end
        k
        E(k) = evidence(end);
        evidence(end)
	end
	toc
	e(m,:)=E;
    e(isnan(e)) = 0;
 	if max(E-E0)<log(sqrt(10))
 		break
 	end
% 	plot(E/E0), drawnow
    BF(m) = max(E-E0);
    if BF > log(sqrt(10))
        [BF(m),k]=max(E-E0);
        J(m)=k;
    end
end
