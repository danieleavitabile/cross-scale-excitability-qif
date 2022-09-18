function [E,k] = ConnectivityMatrix(n,Delta0,kBar)

  Deltak = Delta0*sqrt(kBar);
  k = (kBar+Deltak.*tan((pi/2).*(2.*[1:n]-n-1)./(n+1)))';

  k = k.*( 0<= k & k <=2*kBar);

  k = floor(k);

  inzk = find(k ~= 0);
  knz = k(inzk);

  rows = zeros(sum(knz),1);
  cols = zeros(size(rows));
  vals = ones(size(rows));

  for b = 1:length(knz)

    l = 1+sum(knz(1:b-1));
    u = sum(knz(1:b));
    id = [l:u]';

    iCols = randperm(n,length(id))';
    % pause
    rows(id) = inzk(b);
    cols(id) = randperm(n,length(id))';

  end

  E = sparse(rows,cols,vals,n,n);

end
