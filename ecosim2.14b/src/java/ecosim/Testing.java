package ecosim;

import java.util.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.Date;

public class Testing{

    public static void main(String[] args) {
      
      String newick_file;
      String fasta_file;
      Float cutoff;
      String fasta_out;
      String log_name;
      
      //List<String> settings = import_settings();
      List<String> settings = new ArrayList<String>();
      System.out.println(settings);
      
      if (settings.size() == 5){
        try{
          newick_file = settings.get(0);
          fasta_file = settings.get(1);
          cutoff = Float.parseFloat(settings.get(2));
          fasta_out = settings.get(3);
          log_name = settings.get(4);
        } catch(ArrayIndexOutOfBoundsException e){
          ArrayList<String> logSettings = new ArrayList<String>();
          logSettings.add("Error:Could not find all required settings in 'seq_clean_settings.txt'");
          logSettings.add("-- NEED: newick file, fasta file, cutoff, new fasta file " +
                   "name, new log name. Each should be entered on a new line " +
                   "in the specified order.");
          //logging(logSettings,"SEQ_CLEAN_ERROR_LOG.txt");
        }
      } else{
        ArrayList<String> logFile = new ArrayList<String>();
        logFile.add("Error:Could not find 'seq_clean_settings.txt'");
        //logging(logFile, "SEQ_CLEAN_ERROR_LOG.txt");
        return;
      }
      
      Fasta fasta = new Fasta();

      

      
      NewickTree tree = new NewickTree();//"MLtree from23 Oct with outgroup (fromAligned 31Aug).nwk");
      //NewickTreeNode primary_node = find_node_with_taxon_label("T7C8", tree);

    }

//    private static NewickTreeNode find_node_with_taxon_label(String strain, NewickTree tree){
//        NewickTreeNode nodeFound = new NewickTreeNode();
//        List<NewickTreeNode> descendants = tree.getDescendants();
//        for (int i = 0; i < descendants.size(); i++){
//            NewickTreeNode descendant = descendants.get(i);
//            if (descendant.getName().equals(strain)) {
//                 nodeFound = descendant;
//                 break;
//            }
//        }
//        return nodeFound;
//    }
//    
//    private static NewickTreeNode sister_nodes(NewickTreeNode primary_node){
//      NewickTreeNode sister_found = new NewickTreeNode(); 
//      NewickTreeNode parent_node = primary_node.getParent();
//      List<NewickTreeNode> children = parent_node.getDescendants();
//      for (int i = 0; i < children.size(); i++){
//        NewickTreeNode child = children.get(i);
//        if (!child.getName().equals(primary_node.getName())){
//          sister_found = child;
//          break;
//        }
//      }
//      return sister_found;
//    }
    
    
//    private static List<List<String>> to_tuples(List<String> lst){
//      List<List<String>> outer = new ArrayList<ArrayList<String>>();
//      for (int i = 0; 0 < lst.size() && lst.size() < 2; i++){
//        List<String> inner = new ArrayList<String>();
//        inner.add(lst.get(i));
//        inner.add(lst.get(i+1));
//        outer.add(inner);
//      }
//      return outer;
//    }
//    
//    private static ArrayList<ArrayList<String>> down_search(NewickTreeNode node, Integer dist, ArrayList<ArrayList<String>> lst){
//      if(node.isLeafNode()){
//        //ArrayList<ArrayList> lst = new ArrayList<ArrayList>();
//        ArrayList<String> pair = new ArrayList<String>();
//        pair.add(String.valueOf(dist));
//        pair.add(node.getName());  
//        lst.add(pair);
//        return lst;
//        //return lst.add(pair)
//      }
//      else {
//        dist = dist + 1;
//        ArrayList<NewickTreeNode> children = node.getDescendants();
//        NewickTreeNode child1 = children.get(0);
//        NewickTreeNode child2 = children.get(1);
//        lst = down_search(child1, dist, lst);
//        lst = down_search(child2, dist, lst);
//        // may or may not have to flatten
//        Collections.sort(lst, new Comparator<ArrayList<String>> () {
//          public int compare(ArrayList<String> a, ArrayList<String> b) {
//            return a.get(0).compareTo(b.get(0));
//          }
//        });
//        return lst;
//        //return sorted(to_tuples(flatten([down_search(node.child_nodes()[0],dist),
//       //down_search(node.child_nodes[1],dist)])),key=itemgetter(1))
//      }
//    }
      
//    private static ArrayList<ArrayList<String>> closest_relative(NewickTree tree, String strain){
//     // NewickTreeNode primary_node = find_node_with_taxon_label(strain, tree);
//      NewickTreeNode primary_node = new NewickTreeNode();
//      NewickTreeNode sister_node = new NewickTreeNode();
//      if((primary_node.getName()).equals("")){
//        System.out.println("PRIMARY node " + strain + "could not be found in the tree");
//        return null;
//      }
//      if (primary_node.isLeafNode()){
//        //NewickTreeNode sister_node = sister_nodes(primary_node);
//        ArrayList<ArrayList<String>> lst = new ArrayList<ArrayList<String>>();
//        return down_search(sister_node,0,lst);
//      }
//      else{
//        System.out.println("PRIMARY node " + primary_node + " is NOT a leaf!");
//        return null;
//      }
//    }
//    
//    private static List<String> import_settings(){
//      
//      BufferedReader br = null;
//      List<String> settings = new ArrayList<String>();
//      
//      try{
//        String sCurrentLine;
//        
//        br = new BufferedReader(new FileReader("seq_clean_settings.txt"));
//        
//        while ((sCurrentLine = br.readLine()) != null){
//          settings.add(sCurrentLine);
//        }
//        return settings;
//      } catch (IOException x){
//        return null;
//      }
//    }
    
//    private static void logging(List<String> logs, String log_name){
//      Date date = new Date();
//      String string = "Date: " + date.toString() + "\n" + "\n";
//      for (int i = 0; i < logs.size(); i++){
//        string = string + logs.get(i) + "\n" + "\n";
//      }
//      try{
//        File file = new File(log_name);
//        if (!file.exists()){
//          file.createNewFile();
//        }
//        
//        FileWriter fw = new FileWriter(file.getAbsoluteFile());
//            BufferedWriter bw = new BufferedWriter(fw);
//            bw.write(string);
//            bw.close();
//
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//      }
    
    
}

