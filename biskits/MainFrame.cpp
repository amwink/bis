#include "MainFrame.h"
#include <wx/aboutdlg.h>

MainFrame::MainFrame(wxWindow* parent)
    : MainFrameBaseClass(parent)
{
}

MainFrame::~MainFrame()
{
}

void MainFrame::OnExit(wxCommandEvent& event)
{
    wxUnusedVar(event);
    Close();
}

void MainFrame::OnAbout(wxCommandEvent& event)
{
    wxUnusedVar(event);
    wxAboutDialogInfo info;
    info.SetCopyright(_("Alle Meije Wink"));
    info.SetLicence(_("GPL v2 or later"));
    info.SetDescription(_("BISkits: GUI for the BIS toolset"));
    ::wxAboutBox(info);
}
