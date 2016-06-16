%%
aviobj = avifile('test.avi')
%%
aviobj = close(aviobj);
%%
m = aviread('test.avi');
s = cat(4,m.cdata);
d = 0.1;

for k = 1:size(s,4),
    [a c] = rgb2ind(s(:,:,:,k),256); 
    if k == 1,
        imwrite(a,c,'test.gif','gif','DelayTime',d);
    else        
        imwrite(a,c,'test.gif','gif','writemode','append','DelayTime',d);
    end;
end;
