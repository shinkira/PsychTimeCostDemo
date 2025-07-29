function tg = sample_tguess(obj,tl,th,rand_stream)

if nargin<3
    rand_stream = [];
end

if isempty(rand_stream)
    r = rand(size(th));
else
    r = rand(rand_stream,size(th));
end

tg = tl + (th-tl).*r;


