clear;
%%%%%%%%%%%%%%%%  parameters %%%%%%%%%%%%%%%%%%%%%%
Ne=20056;%31438;%33274;
Nh=43960;%45112;%52456;
c_light=3e8;
k=1.0;
tag='E';
%%%%%%%%%%%%%%% load files  %%%%%%%%%%%%%%%%
load Se_in_in.txt;
Se_in_in=sparse(Se_in_in(:,1)+1,Se_in_in(:,2)+1,Se_in_in(:,3),Nh,Ne);
if tag=='E'
    load Se_b_in.txt;
    Se_b_in=sparse(Se_b_in(:,1)+1,Se_b_in(:,2)+1,Se_b_in(:,3),Nh,Ne);
end
load Sh_in_in.txt;
Sh_in_in=sparse(Sh_in_in(:,1)+1,Sh_in_in(:,2)+1,Sh_in_in(:,3),Ne,Nh);
if tag=='H'
    load Sh_b_in.txt;
    Sh_b_in=sparse(Sh_b_in(:,1)+1,Sh_b_in(:,2)+1,Sh_b_in(:,3),Ne,Nh);
end
load e_in.txt;
e_in=e_in(:,1)+e_in(:,2)*1.0i;
load e_b.txt;
e_b=e_b(:,1)+e_b(:,2)*1.0i;
load h_in.txt;
h_in=h_in(:,1)+h_in(:,2)*1.0i;
load h_b.txt;
h_b=h_b(:,1)+h_b(:,2)*1.0i;
load curlh_in.txt;
curlh_in=curlh_in(:,1)+curlh_in(:,2)*1.0i;
load curle_in.txt;
curle_in=curle_in(:,1)+curle_in(:,2)*1.0i;
load curle_b.txt;
curle_b=curle_b(:,1)+curle_b(:,2)*1.0i;
%%%%%%%%%%%%%%  real code   %%%%%%%%%%%%%%%%%%
if tag=='E'
    Sev=(Se_in_in+Se_b_in)*(e_in+e_b);
    Shv=(Sh_in_in)*(h_in+h_b);
    B=-Sh_in_in*Se_b_in*e_b;%-Sh_b_in*curle_b;
end
if tag=='H'
    Sev=(Se_in_in)*(e_in+e_b);
    Shv=(Sh_in_in+Sh_b_in)*(h_in+h_b);
    B=-Sh_b_in*curle_b;%-Sh_in_in*Se_b_in*e_b;
end

Sev_0=curle_in;
Shv_0=curlh_in;
A=Sh_in_in*Se_in_in-k*k*speye(Ne);


x=A\B;
%sprintf("Error_Se=%f, Error_Sh=%f",norm(Sev-Sev_0)/norm(Sev_0),norm(Shv-Shv_0)/norm(Shv_0))
sprintf("Error_Se=%f, Error_Sh=%f, Error_S=%f",norm(Sev-Sev_0)/norm(Sev_0),norm(Shv-Shv_0)/norm(Shv_0),norm(x-e_in)/norm(e_in))