/* eslint-disable react/self-closing-comp */
import { MemoryRouter as Router, Routes, Route, Link } from 'react-router-dom';
import { useState } from 'react';
import './App.css';

import ImageDisplay from './components/ImageDisplay';
import ElectrodeIPGSelection from './components/ElectrodeIPGSelection';
import TripleToggle from './components/TripleToggle';
import { ReactComponent as MySVG } from './components/electrode_models/images/IPG.svg';
import StimulationParameters from './components/StimulationParameters';
import ElectrodeSide from './components/ElectrodeSide';
import TabbedElectrodeIPGSelection from './components/TabbedElectrodeIPGSelection';
import TabbedElectrodeIPGSelectionTest from './components/TabbedElectrodeIPGSelectionTest';
import MatToJsonConverter from './components/MatToJsonConverter';
import JsonLoaderComponent from './components/extractPolValues';
import JSONDataExtractor from './components/JSONDataExtractor';
import NewBostonCartesiaTest from './components/NewBostonCartesiaTest';
import TripleToggleTest from './components/TripleToggleTest';
import Navbar from './components/Navbar';
import StimulationSettings from './components/StimulationSettings';
import PercentageAmplitudeToggle from './components/PercentageAmplitudeToggle';
import AssistedToggle from './components/AssistedToggle';
import LeadDbsImage from './logo512Padding-300x212.png';
import MAToggleSwitch from './components/MAToggleSwitch';
import ExportData from './components/ExportData';
import React from 'react';
import AssistedButtons from './components/AssistedButtons';
import 'bootstrap/dist/css/bootstrap.min.css';

function Hello() {
  return (
    <div>
      <h1>Stimulation Controller</h1>
      <div className="stimulation-parameters">
        {/* <StimulationParameters /> */}
      </div>
      <div className="stimulation-parameters">
        {/* <StimulationParameters /> */}
      </div>
      <div className="Hello">
        {/* <ElectrodeIPGSelection /> */}
        {/* <JsonLoaderComponent /> */}
        {/* <JSONDataExtractor /> */}
        {/* <TripleToggleTest /> */}
        <TabbedElectrodeIPGSelectionTest />
        {/* <NewBostonCartesiaTest /> */}
        {/* <ElectrodeSide /> */}
        {/* <TripleToggle /> */}
        {/* <StimulationParameters /> */}
        {/* <ElectrodeIPGSelection />
        <ElectrodeModel /> */}
      </div>
    </div>
  );
}

// const RedirectToNewRoute = () => {
//   return <Navigate to="/new-route" />;
// };

export default function App() {
  const [IPG, setIPG] = useState('');
  const [leftElectrode, setLeftElectrode] = useState('');
  const [rightElectrode, setRightElectrode] = useState('');
  // const [key, setKey] = useState('1');
  const [allQuantities, setAllQuantities] = useState({});
  const [allSelectedValues, setAllSelectedValues] = useState({});
  const [allTotalAmplitudes, setAllTotalAmplitudes] = useState({});
  const [allStimulationParameters, setAllStimulationParameters] = useState({});
  const [visModel, setVisModel] = useState('6');
  const [sessionTitle, setSessionTitle] = useState('');
  // const [outputIPG, setOutputIPG] = useState('');

  return (
    <Router>
      <div className="Navbar">
        <Navbar />
        {/* <img src="./logo512Padding-300x212.png" alt="leadDBS" /> */}
      </div>
      <Routes>
        {/* <Route path="/" element={<Hello />} /> */}
        {/* <Route
          path="/"
          element={<Navigate to ="/new-route" />}
        /> */}
        <Route
          path="/"
          element={
            <div>
              <img src={LeadDbsImage} alt="Description of your image" />
              <div></div>
              <Link to="/stimulation-settings">
                <button className="button">Get Started</button>
              </Link>
            </div>
          }
        />
        <Route path="/testing" element={<AssistedButtons />} />
        <Route
          path="/stimulation-settings"
          element={
            <div>
              <StimulationSettings
                IPG={IPG}
                setIPG={setIPG}
                leftElectrode={leftElectrode}
                setLeftElectrode={setLeftElectrode}
                rightElectrode={rightElectrode}
                setRightElectrode={setRightElectrode}
                allQuantities={allQuantities}
                setAllQuantities={setAllQuantities}
                allSelectedValues={allSelectedValues}
                setAllSelectedValues={setAllSelectedValues}
                allTotalAmplitudes={allTotalAmplitudes}
                setAllTotalAmplitudes={setAllTotalAmplitudes}
              />
              <Link to="/tabbed-selection">
                <button className="button">Next</button>
              </Link>
            </div>
          }
        />
        <Route
          path="tabbed-selection"
          element={
            <TabbedElectrodeIPGSelectionTest
              selectedElectrodeLeft={leftElectrode}
              selectedElectrodeRight={rightElectrode}
              IPG={IPG}
              // key={key}
              // setKey={setKey}
              allQuantities={allQuantities}
              setAllQuantities={setAllQuantities}
              allSelectedValues={allSelectedValues}
              setAllSelectedValues={setAllSelectedValues}
              allTotalAmplitudes={allTotalAmplitudes}
              setAllTotalAmplitudes={setAllTotalAmplitudes}
              allStimulationParameters={allStimulationParameters}
              setAllStimulationParameters={setAllStimulationParameters}
              visModel={visModel}
              setVisModel={setVisModel}
              sessionTitle={sessionTitle}
              setSessionTitle={setSessionTitle}
              // outputIPG={outputIPG}
              // setOutputIPG={setOutputIPG}
            />
          }
        />
        <Route
          path="end-session"
          element={
            <ExportData
              allQuantities={allQuantities}
              allSelectedValues={allSelectedValues}
            />
          }
        />
      </Routes>
    </Router>
  );
}
