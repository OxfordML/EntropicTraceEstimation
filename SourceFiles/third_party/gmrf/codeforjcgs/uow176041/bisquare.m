%%%%%%    the bisquare function as a function of distance   %%%%%%%

function R=bisquare(t,r)

% t = distance
% r = range

if numel(r)==1 || isequal(size(r),size(t))
	r_all=r;
else
	if length(vec(r))==size(t,2),
		r_all=repmat(vec(r)',length(t),1);
	else
		disp('ERROR: Bisquare range has wrong dim.')
	end
end

R=(1-(t./r_all).^2).^2.*(t<r_all);
