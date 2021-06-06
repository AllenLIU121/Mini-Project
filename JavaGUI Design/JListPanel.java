package project;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;


public class JListPanel extends JPanel{
	
	private JList list;
	private Vector Items;
	
	public static final int SINGLE_SELECTION=0;
	public static final int SINGLE_INTERVAL_SELECTION=1;
	public static final int MULTIPLE_INTERVAL_SELECTION=2;
	
	public JListPanel(int rows, int height, int width) {
		this(rows, height, width, MULTIPLE_INTERVAL_SELECTION);
	}
	
	public JListPanel(int rows, int height, int width, int choice) {
		list = new JList();
		Items = new Vector();
		
		list.setVisibleRowCount(rows);
		list.setFixedCellHeight(height);
		list.setFixedCellWidth(width);
		if(choice==SINGLE_INTERVAL_SELECTION) {
			list.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		}
		else if(choice==MULTIPLE_INTERVAL_SELECTION) {
			list.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		}
		else
			list.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

		setLayout(new FlowLayout(FlowLayout.CENTER));
		
		add(new JScrollPane(list));
	}
	
	public void addItems(Vector items) {
		Items = items;
		list.setListData(items);
		repaint();
	}	
	
	public void addItems(Object []items) {
		for(int i=0;i<items.length;i++) Items.add(items[i]);
		list.setListData(Items);
		repaint();
	}
	
	public void addItem(Object obj) {
		Items.add(obj);
		list.setListData(Items);
		repaint();
	}
	
	public void removeItem(Object obj) {
		Items.remove(obj);
		list.setListData(Items);
		repaint();
	}
	
	public void removeItems(Object []obj) {
		for(int i=0;i<obj.length;i++)
			removeItem(obj[i]);
	}
	
	public Object[] getSelectedElements() {
		return list.getSelectedValuesList().toArray();
	}
	
	public Object [] getAllElements() {
		Object [] obj = new Object[Items.size()];
		
		for(int i=0; i<Items.size();i++)
			obj[i] = Items.get(i);
		
		return obj;
	}
	
	public void removeAllElements() {
		Items.removeAllElements();
		list.setListData(Items);
		repaint();
	}
	
	public int getElementsCount() {
		return Items.size();
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		JFrame Fr = new JFrame();
		
		JListPanel LP = new JListPanel(5, 20, 100);
		
		for(int i=0;i<20;i++)
			LP.addItem(i);
		
		JButton b = new JButton("Get selected elements");
		b.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				// TODO Auto-generated method stub
				Object []obj = LP.getSelectedElements();
				for(int i=0;i<obj.length;i++)
					System.out.println(obj[i]);
			}
		});	
		
		Fr.getContentPane().add(LP, BorderLayout.CENTER);
		Fr.getContentPane().add(b, BorderLayout.SOUTH);
		
		Fr.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		Fr.pack();
		Fr.setVisible(true);
	}

}
