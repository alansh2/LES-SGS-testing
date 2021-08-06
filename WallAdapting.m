function ev = WallAdapting(Delta,S,g);
mkdir temp
Sid = {'S11','S12','S13';
       'S12','S22','S23';
       'S13','S23','S33'};
idx = {'i11','i12','i13';
       'i21','i22','i23';
       'i31','i32','i33'};
sz = size(g,'dUdx');
S2 = S.S11.^2 + S.S22.^2 + S.S33.^2 + 2*(S.S12.^2 + S.S13.^2 + S.S23.^2);
O = matfile('./temp/Omega.mat','Writable',true);
O.i11 = zeros(sz); O.i22 = O.i11; O.i33 = O.i22;
O.i12 = 0.5*(g.dUdy - g.dVdx); O.i21 = -O.i12;
O.i13 = 0.5*(g.dUdz - g.dWdx); O.i31 = -O.i13;
O.i23 = 0.5*(g.dVdz - g.dWdy); O.i32 = -O.i23;
O2 = 2*(O.i12.^2 + O.i13.^2 + O.i23.^2);
SS = matfile('./temp/SikSkj.mat','Writable',true);
OO = matfile('./temp/OjlOli.mat','Writable',true);
for i = 1:3
    for j = 1:3
        SS.(idx{i,j}) = S.(Sid{i,1}).*S.(Sid{1,j}) + S.(Sid{i,2}).*S.(Sid{2,j}) + S.(Sid{i,3}).*S.(Sid{3,j});
        OO.(idx{i,j}) = O.(idx{i,1}).*O.(idx{1,j}) + O.(idx{i,2}).*O.(idx{2,j}) + O.(idx{i,3}).*O.(idx{3,j});
    end
end
IV = zeros(sz);
for i = 1:3
    IV = IV + SS.(idx{i,1}).*OO.(idx{1,i}) + SS.(idx{i,2}).*OO.(idx{2,i}) + SS.(idx{i,3}).*OO.(idx{3,i});
end
Sd2 = 1/6*(S2.^2 + O2.^2) + 2/3*S2.*O2 + 2*IV;
cw = 0.56;
ev = (cw*Delta)^2*(Sd2.^1.5)./(S2.^2.5 + Sd2.^1.25);
delete('./temp/*.mat');
