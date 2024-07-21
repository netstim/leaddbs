import express from 'express';
import bodyParser from 'body-parser';
import cors from 'cors'; // Import cors middleware

const app = express();
const port = process.env.PORT || 3001;

app.use(bodyParser.json());
app.use(cors()); // Enable CORS for all routes

let receivedData: any;

// Define your API routes
app.get('/api/data', (req, res) => {
  // Handle GET request for /api/data
  // res.json({ message: 'Data from server' });
  if (receivedData) {
    // res.json({ message: 'Data sent' });
    res.json(receivedData);
  } else {
    res.json({ message: 'Data not sent' });
  }
});

// Handle POST request to /api/data
app.post('/api/data', (req, res) => {
  // Handle POST request for /api/data
  receivedData = req.body;
  console.log(req.body); // Log the data received from the client
  res.json({ message: 'Data received by server' });
});

// Start the server
app.listen(port, () => {
  console.log(`Server running at http://localhost:${port}`);
});
