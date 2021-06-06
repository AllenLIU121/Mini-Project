package project;

import javax.swing.*;
import java.awt.event.*;

public class SimpleMenu extends JMenuBar{
	public SimpleMenu(String MB[][], ActionListener al) {
		for(int i=0; i<MB.length; i++) createMenu(MB[i], al);
	}
	
	private void createMenu(String M[], ActionListener al) {
		JMenu jm = new JMenu(M[0]);
		for (int i=1; i<M.length; i++) {
			if (M[i].equals("-")) jm.addSeparator();
			else {
				JMenuItem mi = new JMenuItem(M[i],new ImageIcon("images/"+M[i]+".jpg"));
				mi.addActionListener(al);
				mi.setMnemonic(M[i].charAt(0));
				jm.add(mi);
			}
		}
		add(jm);
	}
}