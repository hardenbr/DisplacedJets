// Joshua Hardenbrook 18-11-2015
// 1. encapsulates the event level selection, as well as the jet variable selection
// 2. parses the variables / cuts from a JSON


class jetSelector {

 public:

  // parse the input selection 
  void parseSelectionJSON();

  // add a requirement on the events used in the selection
  void addEventReq();
  // add a requirement on the jets to be selected
  void addJetReq();

  // prints a summary of all the current requirements on events and jets
  void printSelectionSummary();
  
  // check if a jet is selectd
  void isSelected(const jetCandidate&);

  // indices for the number of requirements in the event
  int nJetRequirements;
  int nEventRequirements;    

}
