function [th] = sample_th(Th_vert)
%SAMPLE_TH sample a point th in the polytope formed by Th_vert
% this is not the uniform distribution, it is just a simple sampler.
c = rand(size(Th_vert, 2), 1);

cn = abs(c)/norm(c,1);

th = Th_vert*cn;
end

