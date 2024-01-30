function ea_conda_setproxy

proxy = inputdlg('Set proxy for conda and pip', '', 1, {'http://'});
proxy = proxy{1};
ea_cprintf('*Comments', 'Set proxy for pip...\n');
ea_conda.run(['pip config set global.proxy ' proxy]);
ea_cprintf('*Comments', 'Set proxy for conda...\n');
ea_conda.run(['conda config --set proxy_servers.http ' proxy]);
ea_conda.run(['conda config --set proxy_servers.https ' proxy]);
