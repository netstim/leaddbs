/* eslint global-require: off, no-console: off, promise/always-return: off */

/**
 * This module executes inside of electron's main process. You can start
 * electron renderer process from here and communicate with the other processes
 * through IPC.
 *
 * When running `npm run build` or `npm run build:main`, this file is compiled to
 * `./src/main.js` using webpack. This gives us some performance wins.
 */
import path from 'path';
import { app, BrowserWindow, shell, ipcMain } from 'electron';
import { autoUpdater } from 'electron-updater';
import log from 'electron-log';
import * as childProcess from 'child_process';
import MenuBuilder from './menu';
import { resolveHtmlPath } from './util';

class AppUpdater {
  constructor() {
    log.transports.file.level = 'info';
    autoUpdater.logger = log;
    autoUpdater.checkForUpdatesAndNotify();
  }
}

let mainWindow: BrowserWindow | null = null;

const startServer = () => {
  // Start the Express server in a child process
  const serverProcess = childProcess.spawn('node', ['dist/server.js'], {
    cwd: path.join(__dirname, '../'), // Adjust the path as needed
    stdio: 'inherit',
  });

  serverProcess.on('error', (err) => {
    console.error('Failed to start server:', err);
  });

  serverProcess.on('exit', (code, signal) => {
    console.log('Server process exited with code:', code);
  });

  return serverProcess;
};

app
  .whenReady()
  .then(() => {
    // Start Express server
    const serverProcess = startServer();

    // Create window and other initialization code...
  })
  .catch(console.error);

ipcMain.on('ipc-example', async (event, arg) => {
  const msgTemplate = (pingPong: string) => `IPC test: ${pingPong}`;
  console.log(msgTemplate(arg));
  const fs = require('fs');

  // const f = fs.readFileSync('/Users/savirmadan/Documents/GitHub/leaddbs/tempData.json');
  const f = fs.readFileSync(
    '/Users/savirmadan/Development/lead-dbs-programmer/tempData.json',
  );
  console.log(f);
  // const separatedLine = f.split('\\\\');
  // console.log(separatedLine);
  // const k = fs.readFileSync(f);
  // console.log(k);
  event.reply('ipc-example', msgTemplate(`pong: ${f}`));
});

ipcMain.on('import-file', async (event, arg) => {
  const msgTemplate = (pingPong: string) => `${pingPong}`;
  const fs = require('fs');
  const currentDirectory = app.getAppPath();
  const directories = currentDirectory.split('/');

  // Initialize a variable to store the result
  let result = '';

  // Loop through the directories
  for (const dir of directories) {
    // Append each directory to the result
    result += `${dir}/`;

    // If the directory contains "lead-dbs-programmer", stop the loop
    if (dir === 'lead-dbs-programmer') {
      break;
    }
  }

  const fileName = 'inputData.json';
  const filePath = path.join(result, fileName);
  console.log(filePath);
  const f = fs.readFileSync(filePath);
  const jsonData = JSON.parse(f);
  event.reply('import-file', jsonData);
});

ipcMain.on('open-file', (event, arg) => {
  // const msgTemplate = (pingPong: string) => `IPC test: ${pingPong}`;
  // console.log(msgTemplate(arg));
  const fs = require('fs');

  // const f = fs.readFileSync('/Users/savirmadan/Documents/GitHub/leaddbs/tempData.json');
  const f = fs.readFileSync(arg);
  console.log(event);
  // const separatedLine = f.split('\\\\');
  // console.log(separatedLine);
  // const k = fs.readFileSync(f);
  // console.log(k);
  event.reply('open-file', `pong: ${f}`);
});

// ipcMain.on('close-window', () => {
//   // const currentWindow = BrowserWindow.getFocusedWindow();
//   // if (currentWindow) {
//   //   currentWindow.close();
//   //   event.reply('window-closed', 'Window closed');
//   // }

//   app.quit();
// });

ipcMain.on('close-window', (event, arg) => {
  // setTimeout(() => {
  //   app.quit();
  // }, 5000); // 5000 milliseconds = 5 seconds

  app.quit();

  // const currentWindow = BrowserWindow.getFocusedWindow();
  // if (currentWindow) {
  //   currentWindow.close();
  //   event.reply('window-closed', 'Window closed');
  // }
});

const { dialog } = require('electron');
const fs = require('fs');

ipcMain.on('save-file', (event, data) => {
  // Example of saving data to a file
  // const filePath = app.getPath('downloads') + '/data.txt';
  // const currentDirectory = app.getAppPath();
  // const currentDirectory = '/Users/savirmadan/Development/lead-dbs-programmer';
  const currentDirectory = app.getAppPath();
  const directories = currentDirectory.split('/');

  // Initialize a variable to store the result
  let result = '';

  // Loop through the directories
  for (const dir of directories) {
    // Append each directory to the result
    result += `${dir}/`;

    // If the directory contains "lead-dbs-programmer", stop the loop
    if (dir === 'lead-dbs-programmer') {
      break;
    }
  }
  // console.log(currentDirectory + '/lead-dbs-programmer');
  if (currentDirectory) {
    // Convert data to string format
    const dataString = JSON.stringify(data);
    const fileName = 'data.json';
    const filePath = path.join(result, fileName);
    // const filePath = './dist/main/webpack:/leaddbs-stimcontroller/main.js';
    // Write data to file
    fs.writeFileSync(filePath, dataString);

    // Send a response back to the renderer process
    event.reply('file-saved', filePath);
  }
});

// ipcMain.on('open-file', (event, filePath) => {
//   // Read file data
//   const fs = require('fs');
//   fs.readFile(filePath, 'utf8', (err: any, data: any) => {
//     if (err) {
//       // Handle error
//       console.error(err);
//       event.sender.send('file-data', null);
//       return;
//     }
//     // Send data back to renderer process
//     event.sender.send('file-data', data);
//   });
// });

// const { spawn } = require('child_process');

// ipcMain.on('trigger-matlab-action', (event, data) => {
//   const matlabProcess = spawn('matlab', ['-r', 'disp("hello")']);

//   matlabProcess.stdout.on('data', (data) => {
//       console.log(`MATLAB stdout: ${data}`);
//   });

//   matlabProcess.stderr.on('data', (data) => {
//       console.error(`MATLAB stderr: ${data}`);
//   });
// });

if (process.env.NODE_ENV === 'production') {
  const sourceMapSupport = require('source-map-support');
  sourceMapSupport.install();
}

const isDebug =
  process.env.NODE_ENV === 'development' || process.env.DEBUG_PROD === 'true';

if (isDebug) {
  require('electron-debug')();
}

const installExtensions = async () => {
  const installer = require('electron-devtools-installer');
  const forceDownload = !!process.env.UPGRADE_EXTENSIONS;
  const extensions = ['REACT_DEVELOPER_TOOLS'];

  return installer
    .default(
      extensions.map((name) => installer[name]),
      forceDownload,
    )
    .catch(console.log);
};

const createWindow = async () => {
  if (isDebug) {
    await installExtensions();
  }

  const RESOURCES_PATH = app.isPackaged
    ? path.join(process.resourcesPath, 'assets')
    : path.join(__dirname, '../../assets');

  const getAssetPath = (...paths: string[]): string => {
    return path.join(RESOURCES_PATH, ...paths);
  };

  mainWindow = new BrowserWindow({
    show: false,
    width: 1100,
    height: 1100,
    maxWidth: 1100, // Maximum width of the window
    // maxHeight: 1200, // Maximum height of the window
    minWidth: 1000, // Minimum width of the window
    // minHeight: 1200, // Minimum height of the window
    icon: getAssetPath('icon.png'),
    webPreferences: {
      preload: app.isPackaged
        ? path.join(__dirname, 'preload.js')
        : path.join(__dirname, '../../.erb/dll/preload.js'),
    },
  });

  mainWindow.loadURL(resolveHtmlPath('index.html'));

  mainWindow.on('ready-to-show', () => {
    if (!mainWindow) {
      throw new Error('"mainWindow" is not defined');
    }
    if (process.env.START_MINIMIZED) {
      mainWindow.minimize();
    } else {
      mainWindow.show();
    }
  });

  mainWindow.on('closed', () => {
    mainWindow = null;
  });

  const menuBuilder = new MenuBuilder(mainWindow);
  menuBuilder.buildMenu();

  // Open urls in the user's browser
  mainWindow.webContents.setWindowOpenHandler((edata) => {
    shell.openExternal(edata.url);
    return { action: 'deny' };
  });

  // Remove this if your app does not use auto updates
  // eslint-disable-next-line
  new AppUpdater();
};

/**
 * Add event listeners...
 */

app.on('window-all-closed', () => {
  // Respect the OSX convention of having the application in memory even
  // after all windows have been closed
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app
  .whenReady()
  .then(() => {
    createWindow();
    app.on('activate', () => {
      // On macOS it's common to re-create a window in the app when the
      // dock icon is clicked and there are no other windows open.
      if (mainWindow === null) createWindow();
    });
  })
  .catch(console.log);
