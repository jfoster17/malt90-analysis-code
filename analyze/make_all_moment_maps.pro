pro make_all_moment_maps
lines = ['n2hp','13cs','h41a','ch3cn','hc3n','13c34s','hnc','hc13ccn','hcop','hcn','hnco413','hnco404','c2h','sio','hn13c','h13cop']
for i=0,n_elements(lines)-1 do begin
    directory = lines[i]
    moment_maps_malt,directory
endfor
end
