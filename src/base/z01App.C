#include "z01App.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
z01App::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

z01App::z01App(InputParameters parameters) : MooseApp(parameters)
{
  z01App::registerAll(_factory, _action_factory, _syntax);
}

z01App::~z01App() {}

void
z01App::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<z01App>(f, af, s);
  Registry::registerObjectsTo(f, {"z01App"});
  Registry::registerActionsTo(af, {"z01App"});

  /* register custom execute flags, action syntax, etc. here */
}

void
z01App::registerApps()
{
  registerApp(z01App);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
z01App__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  z01App::registerAll(f, af, s);
}
extern "C" void
z01App__registerApps()
{
  z01App::registerApps();
}
