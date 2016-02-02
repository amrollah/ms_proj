function sendmail(from,to,subject,text)
import java.io.*;
import java.net.*;
try
  u = URL(['mailto:' to]); % Create a mailto: URL
  c = u.openConnection; % Create its URLConnection
  c.setDoInput(0); % Specify no input from it
  c.setDoOutput(1); % Specify we'll do output
  c.connect; % Connect to mail host
  osw = OutputStreamWriter(c.getOutputStream);
  out = PrintWriter(osw);
  user = char(java.lang.System.getProperty('user.name'));
  host = char(InetAddress.getLocalHost.getHostName);
  out.print(['From: "' from '" <' user '@' host '>\n']);
  out.print(['To: ' to '\n']);
  out.print(['Subject: ' subject '\n']);
  out.print('\n');  % blank line to end the list of headers
  if ischar(text) out.print(text); 
  elseif isnumeric(text) out.print(num2str(text));  
  elseif iscell(text)
    bWarned = 0;
    text = text{:};
    for i=1:length(text)
      if ischar(text{i}) out.print(text{i}); 
      elseif isnumeric(text{i}) out.print(num2str(text{i}));
      elseif ~bWarned 
        warning(['I dont know how to email objects of type ' class(text{i})]);
        bWarned=1;
      end;
      out.print('\n'); 
    end;
  else 
    warning(['I dont know how to email objects of type ' class(text)]);
  end;
  out.print('\n'); 
  out.close;
catch
  error('sendmail failed');
end;