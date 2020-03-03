#include "Fortran/UAPFortran.hpp"

#include "UAP/UAPUtilities.hpp"
#include "AML/AMLReader.hpp"
#include "AML/AMLLatticeExpander.hpp"
#include "Translate/BmadParser.hpp"

using namespace std;

// Routines used by aml_parser_cpp

class UAPToBmadLat : public BmadParser {

public:

  void init (UAPNode* expand_root);
  void register_aml_element (UAPNode* parent);
  void aml_element_to_bmad (UAPNode* parent);
  void aml_parameter_to_bmad (string who, UAPNode* expand_node, UAPNode* bmad_param_node);
  void aml_parser_cpp(char* lat_file, uap_node_f_struct* fort_root, int& ok);

};

//--------------------------------------------------------------------------------------
// void aml_parser_cpp_(lat_file, root)
//

extern "C" void aml_parser_cpp_(char* lat_file, uap_node_f_struct* fort_root, int& ok) {
  UAPToBmadLat ubl;
  ubl.aml_parser_cpp (lat_file, fort_root, ok);
  return;
}

//-----------------------------------------------------------------------------

void UAPToBmadLat::aml_parser_cpp (char* lat_file, uap_node_f_struct* fort_root, int& ok) {

  UAPNode* uap_root_node;
  UAPNode* expand_root;

  // Read in the AML file and create an AML representation tree.

  ok = false;

  AMLReader reader;
  try {
    uap_root_node = reader.AMLFileToAMLRep (lat_file);
    if (!uap_root_node) return;
  } catch (UAPFileCouldNotBeOpened err) {
    cerr << err << endl;
    return;
  }

  // Expand the lattice

  try {
    AMLReader reader;
    AMLLatticeExpander LE;
    expand_root = LE.AMLExpandLattice(uap_root_node);
    if (!expand_root) return;
  } catch (UAPException err) {
    cerr << err << endl;
    return;
  }

  // Output expanded lattice.

  ofstream aml_out("expanded_lat.aml");
  aml_out << expand_root->toStringTree();
  cout << "Expanded AML Lattice written to: expanded_lat.aml" << endl;
  aml_out.close();

  // Put in the appropriate Bmad class into each element

  UAPNode* track_node = expand_root->getSubNodesByName("tracking_lattice").front();
  UAPNode* master_node = expand_root->getSubNodesByName("master_list").front();
  UAPNode* control_node = expand_root->getSubNodesByName("control_list").front();

  init (expand_root);

  register_aml_element (track_node);
  register_aml_element (master_node);
  register_aml_element (control_node);

  aml_element_to_bmad (track_node);
  aml_element_to_bmad (master_node);
  aml_element_to_bmad (control_node);

  UAPNode* param_node = expand_root->addChild("param_list");

  aml_parameter_to_bmad ("beam", expand_root, param_node);
  aml_parameter_to_bmad ("lattice", expand_root, param_node);
  aml_parameter_to_bmad ("global", expand_root, param_node);

  // Convert to a fortran tree.

  uap_node_tree_to_f_(expand_root, fort_root);
  ok = true;

}

//-----------------------------------------------------------------------------

void UAPToBmadLat::init (UAPNode* expand_root) {
  aml_rep_to_x_init (expand_root);
  return;
}

//-----------------------------------------------------------------------------

void UAPToBmadLat::register_aml_element (UAPNode* parent) {

  info_out.ix_line = -1;
  info_out.parsing_status = "Converting AML lattice structures to Bmad.";

  // This is necessary so the call to aml_node_to_x_rep (called in 
  // aml_element_to_x_class) can find the bmad_class of any slave element of a controller.

  NodeVec children = parent->getChildren();
  for (NodeVecIter in = children.begin(); in != children.end(); in++) {
    UAPNode* ele_node = *in;
    string bmad_class;
    string ele_name = ele_node->getAttributeString("name");
    aml_ele_to_x_class (ele_node, bmad_class);
    name_to_x_class_map[ele_name] = bmad_class;
  }
}

//-----------------------------------------------------------------------------

void UAPToBmadLat::aml_element_to_bmad (UAPNode* parent) {

  NodeVec children = parent->getChildren();
  for (NodeVecIter in = children.begin(); in != children.end(); in++) {
    UAPNode* ele_node = *in;

    // Find the Bmad class (key) for this element

    string bmad_class;
    aml_ele_to_x_class (ele_node, bmad_class);

    // Translate the element attribute names.
    // We put the translation into a node and at the end attach it to the element node.

    UAPNode* bmad_attrib_node = new UAPNode("bmad_attributes");

    if (bmad_class == "group" || bmad_class == "overlay") {
      aml_node_to_x (ele_node, bmad_attrib_node);
      UAPNode* node = bmad_attrib_node->getChildren().back();
      node = ele_node->addChildCopy(node);
      node->addAttribute("class", node->getName());
      node->setName("bmad_attributes");
      bmad_attrib_node->deleteNode();
    } else {
      bmad_attrib_node->addAttribute("class", bmad_class);    
      aml_element_attributes_to_x (ele_node, ele_node, bmad_attrib_node);
      ele_node->addChildCopy(bmad_attrib_node);
      bmad_attrib_node->deleteNode();
    }

  }

}

//-----------------------------------------------------------------------------

void UAPToBmadLat::aml_parameter_to_bmad (string who, UAPNode* expand_node, 
                                                    UAPNode* bmad_param_node) {

  if (expand_node->getSubNodesByName(who).size() == 0) return;
  UAPNode* aml_param_node = expand_node->getSubNodesByName(who).front();
  aml_param_to_x (aml_param_node, aml_param_node, bmad_param_node);

}
