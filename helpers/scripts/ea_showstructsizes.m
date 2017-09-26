function ea_showstructsizes(S)

fn=fieldnames(S);
for f=1:length(fn)
   A=S.(fn{f});
   
   R=whos('A');
   disp([fn{f},': ',sprintf('%0.1f',R.bytes/1000),' KB.']);    
end