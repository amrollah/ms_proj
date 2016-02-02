function getfrombase(vname)
assignin('caller',vname,evalin('base',vname));
