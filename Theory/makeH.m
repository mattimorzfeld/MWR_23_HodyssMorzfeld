function H = makeH(no,nr,nx)
ny = no*nr;
H = zeros(ny,nx);
for ll=1:no
    h = zeros(no*nr,1);
    h(1+(ll-1)*nr:ll*nr) = 1;
    H(:,ll) = h;
end