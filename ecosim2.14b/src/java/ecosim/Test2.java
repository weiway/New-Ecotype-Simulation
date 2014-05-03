package ecosim;

import java.util.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.IOException;

public class Test2{ 

  public static void main(String[] args) throws IOException, InvalidNewickException{
    
    String newick_file;
    String fasta_file;
    double cutoff;
    String fasta_out;
    String log_name;
      
    // importing settings from a text file called 'seq_clean_settings.txt'
    List<String> settings = import_settings();
    if(settings.size() == 5){
      try{
        newick_file = settings.get(0);
        fasta_file = settings.get(1);
        cutoff = Double.parseDouble(settings.get(2));
        fasta_out = settings.get(3);
        log_name = settings.get(4);
      } catch(ArrayIndexOutOfBoundsException e){
        ArrayList<String> logSettings = new ArrayList<String>();
        logSettings.add("Error:Could not find all required settings in 'seq_clean_settings.txt'");
        logSettings.add("-- NEED: newick file, fasta file, cutoff, new fasta file " +
                        "name, new log name. Each should be entered on a new line " +
                        "in the specified order.");
        logging(logSettings,"SEQ_CLEAN_ERROR_LOG.txt");
        return;
      }} else{
        ArrayList<String> logFile = new ArrayList<String>();
        logFile.add("Error:Could not find 'seq_clean_settings.txt'");
        logging(logFile, "SEQ_CLEAN_ERROR_LOG.txt");
        return;
      }
      // load in fasta file and get strain:sequences in a hashmap
      Fasta fasta = new Fasta(fasta_file);
      HashMap<String,String> fasta_map = fasta.getSequenceHash();
      HashMap<String,String> fasta_map_copy = (HashMap<String,String>)fasta_map.clone();
      HashMap<String,String> new_fasta_map = replacer(newick_file,fasta_map_copy,cutoff);
      update_fasta(fasta,new_fasta_map);
      fasta.save(fasta_out);

  }
  /**
   * Returns the NewickTreeNode with the specified taxon label in the given tree.
   * 
   * @param strain Name of strain being searched for.
   * @param tree Newick formatted tree to search down.
   * @return NewickTreeNode containing the wanted strain from the Newick formatted tree.
   */
  private static NewickTreeNode find_node_with_taxon_label(String strain, NewickTree tree){
    NewickTreeNode nodeFound = new NewickTreeNode();
    List<NewickTreeNode> descendants = tree.getDescendants();
    for (int i = 0; i < descendants.size(); i++){
      NewickTreeNode descendant = descendants.get(i);
      if (descendant.getName().equals(strain)) {
        nodeFound = descendant;
        break;
      }
    }
        return nodeFound;
  }
  
  private static List<NewickTreeNode> allDescendants(NewickTreeNode node){
    ArrayList<NewickTreeNode> descendants = new ArrayList<NewickTreeNode>();
    ArrayList<NewickTreeNode> children = node.getChildren();
    if (children.size() > 0) {
      for (int i = 0; i < children.size(); i ++) {
        NewickTreeNode child = children.get(i);
        descendants.add(child);
      }
    }
    return descendants;
  }
    
    
  /**
   * Returns the sister of the NewickTreeNode given.
   * 
   * @param primary_node NewickTreeNode to find the sister of.
   * @return NewickTreeNode containing the sister of the node given.
   */
   private static NewickTreeNode sister_nodes(NewickTreeNode primary_node){
      NewickTreeNode sister_found = new NewickTreeNode(); 
      NewickTreeNode parent_node = primary_node.getParent();
      List<NewickTreeNode> children = allDescendants(parent_node);
      for (int i = 0; i < children.size(); i++){
        NewickTreeNode child = children.get(i);
        if (!child.getName().equals(primary_node.getName())){
          sister_found = child;
          break;
        }
      }
      return sister_found;
    }
  
  /**
   * Returns a sorted list of [[name,distance],[name2,distance2]] based on distance.
   * 
   * @param node NewickTreeNode containing the beginning node to search?
   * @param dist int containing the distance from the original node.
   * @param lst List<List<String>> containing [[name,distance]]
   * @return List<List<String>> containing [[name,distance]]
   */
  private static List<List<String>> down_search(NewickTreeNode node, int dist, List<List<String>> lst){
      if(node.isLeafNode()){
        List<String> pair = new ArrayList<String>();
        pair.add(node.getName());  
        pair.add(String.valueOf(dist));
        // removed from list?
        lst.add(pair);
        return lst;
        //return lst.add(pair)
      }
      else {
        dist = dist + 1;
        List<List<String>> down = new ArrayList<List<String>>();
        List<NewickTreeNode> children = allDescendants(node);
        
        for (int i = 0; i < children.size(); i++){
          NewickTreeNode child = children.get(i);
          down_search(child,dist,lst);
        }
        
        // may or may not have to flatten
        Collections.sort(lst, new Comparator<List<String>> () {
          public int compare(List<String> a, List<String> b) {
            return a.get(1).compareTo(b.get(1));
          }
        });
        return lst;
      }
  }
  
  /**
   * Returns a sorted list of [[name,distance],[name2,distance2]] based on distance, 
   * indicating the closest relatives.
   * 
   * @param tree Newick formatted tree to search.
   * @param strain Name of the strain's relatives to find.
   * @return List<List<String>> containing [[name,distance]] in order of closest relatives. 
   */
  private static List<List<String>> closest_relative(NewickTree tree, String strain){
      NewickTreeNode primary_node = find_node_with_taxon_label(strain, tree);
      if((primary_node.getName()).equals("")){
        System.out.println("PRIMARY node " + strain + "could not be found in the tree");
        return null;
      }
      if (primary_node.isLeafNode()){
        NewickTreeNode sister_node = sister_nodes(primary_node);
        List<List<String>> lst = new ArrayList<List<String>>();
        return down_search(sister_node,0,lst);
      }
      else{
        System.out.println("PRIMARY node " + primary_node + " is NOT a leaf!");
        return null;
      }
    }
  
  private static HashMap<NewickTreeNode,List<List<String>>> get_parents_parent_relatives(NewickTreeNode strain_node){
    HashMap<NewickTreeNode,List<List<String>>> rel = new HashMap<NewickTreeNode,List<List<String>>>();
    NewickTreeNode parent2 = strain_node.getParent().getParent();
    List<List<String>> lst = new ArrayList<List<String>>();
    rel.put(parent2, down_search(parent2,0,lst));
    return rel;
  }
  
  /** Imports settings from file. Expects a newick file, a fasta file,
   *  a name for the new fasta file, and a name for the log file
   *  THE SETTINGS FILE MUST BE CALLED 'seq_clean_settings.txt' and must be in
   *  the current directory! If not found, will write to log file that it can't find
   *  settings in CWD.
   *  
   *  @return List<String> containing the 5 settings from the file.
   */
    private static List<String> import_settings(){
      
      BufferedReader br = null;
      List<String> settings = new ArrayList<String>();
      
      try{
        String sCurrentLine;
        
        br = new BufferedReader(new FileReader("seq_clean_settings.txt"));
        
        while ((sCurrentLine = br.readLine()) != null){
          settings.add(sCurrentLine);
        }
        return settings;
      } catch (IOException x){
        return null;
      }
    }
    
      private static void logging(List<String> logs, String log_name){
      Date date = new Date();
      String string = "Date: " + date.toString() + "\n" + "\n";
      for (int i = 0; i < logs.size(); i++){
        string = string + logs.get(i) + "\n" + "\n";
      }
      try{
        File file = new File(log_name);
        if (!file.exists()){
          file.createNewFile();
        }
        
        FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write(string);
            bw.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
      }
      
      private static double percent_gap(HashMap<String,String> fasta_map, int index){
    int num_strains = fasta_map.size();
    int gaps = 0;
    Iterator it = fasta_map.entrySet().iterator();
    while (it.hasNext()){
      Map.Entry pairs = (Map.Entry)it.next();
      String value = (String)pairs.getValue();
      try{
      if (value.charAt(index) == '-'){
          gaps = gaps + 1;
        }} catch(IndexOutOfBoundsException e) {
          gaps = gaps + 1;
        }
    }
    return (double) gaps/num_strains;
  }
  
  private static HashMap<String,String> remove_column(HashMap<String,String> fasta_map, 
                                                      int index){
    HashMap<String,String> new_fasta_map = new HashMap<String,String>();
    Iterator it = fasta_map.entrySet().iterator();
    while (it.hasNext()){
      Map.Entry pairs = (Map.Entry)it.next();
      String key = (String)pairs.getKey();
      String value = (String)pairs.getValue();
      List<String> vlist = new ArrayList<String>(Arrays.asList(value.split("")));
      vlist.remove(0);
      try{
        vlist.remove(index);
      } catch(IndexOutOfBoundsException e) {
        ;
      }
      StringBuilder sb = new StringBuilder();
      for (String str : vlist){
        sb.append(str);
      }
      String new_value = sb.toString();
      new_fasta_map.put(key,new_value);
    }
    return new_fasta_map;
  }
  // initialize level = 0
  @SuppressWarnings("unchecked")
  private static HashMap<Character,List<List<List<String>>>> weighted_seq_choose(HashMap<String,String> strain_fasta_map, 
                                     NewickTreeNode strain_node, 
                                     List<List<List<String>>> rels, 
                                     int index, int level){
    HashMap<Character,Integer> nucs = new HashMap<Character,Integer>();
    char[] nuc_list =  new char[] {'A','T','C','G'};
    List<List<String>> rel_level = rels.get(level); 
    for(List<String> r : rel_level){
      try{
        char nucleotide = (strain_fasta_map.get(r.get(0))).charAt(index);
        if(nucs.containsKey(nucleotide)){
          int value = nucs.get(nucleotide);
          nucs.put(nucleotide,value+1);
        }else{
          nucs.put(nucleotide,1);
        }
      } catch(IndexOutOfBoundsException e) {
        ;
      }
    }
    List<Character> pop = new ArrayList<Character>();
    Iterator it = nucs.entrySet().iterator();
    while (it.hasNext()){
      Map.Entry pairs = (Map.Entry)it.next();
      char key = (Character)pairs.getKey();
      int value = (Integer)pairs.getValue();
      if(Arrays.asList(nuc_list).contains(key)){
        for(int i = 0; i < value; i++){
          pop.add(key);
        }
      }
    }
      if(pop.isEmpty()) {
        List<List<String>> rels_next = new ArrayList<List<String>>();
        NewickTreeNode parent2 = new NewickTreeNode();
        try{
          rels_next = rels.get(level+1); 
          parent2 = strain_node.getParent();
        }catch(IndexOutOfBoundsException e){
          HashMap<NewickTreeNode,List<List<String>>> rels_next_map = 
            get_parents_parent_relatives(strain_node);
          Iterator it2 = rels_next_map.entrySet().iterator();
          while (it2.hasNext()){
            Map.Entry pairs = (Map.Entry)it2.next();
            parent2 = (NewickTreeNode)pairs.getKey();
            rels_next = (List<List<String>>)pairs.getValue(); 
            rels.add(rels_next);
        }
      }
        return weighted_seq_choose(strain_fasta_map,parent2,rels,index,level+1);
      }
      Random rand = new Random();
      int rand_num = rand.nextInt(pop.size());
      char choice = pop.get(rand_num);
      HashMap<Character,List<List<List<String>>>> char_rels = new HashMap<Character,List<List<List<String>>>>();
      char_rels.put(choice,rels);
      return char_rels;
  }
  
  @SuppressWarnings("unchecked")
  private static HashMap<String,String> replacer(String newick_file, 
                                                 HashMap<String,String> fasta_map, 
                                                 double cutoff) throws IOException, InvalidNewickException{
    String stree = new Scanner(new File(newick_file)).useDelimiter("\\z").next();
    NewickTree tree = new NewickTree(stree);
    HashMap<String,String> orig_fasta_map = (HashMap<String,String>)fasta_map.clone();
    HashMap<String,String> new_fasta_map = new HashMap<String,String>();
    char[] nuc_list =  new char[] {'A','T','C','G'};
    int s_num = 0;
    int c_removals = 0;
    int g_replacements = 0;
    int base_replacements = 0;
    String closest = "";
    System.out.println(orig_fasta_map.keySet());
    
    for(String strain : fasta_map.keySet()){
      int index = 0;
      List<List<List<String>>> rels = new ArrayList<List<List<String>>>();
      rels.add(closest_relative(tree,strain));
      String sequence_copy = fasta_map.get(strain);
      NewickTreeNode strain_node = new NewickTreeNode();
      
      for(int i = 0; i < sequence_copy.length(); i++){
        char nuc = sequence_copy.charAt(i);
        if(nuc == '-'){
          double p_gaps = percent_gap(fasta_map,index);
          if(p_gaps <= cutoff){
            if(strain_node == null){
              strain_node = find_node_with_taxon_label(strain,tree);
            }
            
            HashMap<Character,List<List<List<String>>>> close_rel = 
              weighted_seq_choose(orig_fasta_map,strain_node,rels,index,0);
            Iterator it = close_rel.entrySet().iterator();
            while (it.hasNext()){
              Map.Entry pairs = (Map.Entry)it.next();
              closest = (String)pairs.getKey();
              rels = (List<List<List<String>>>)pairs.getValue();
            }
            if(closest != ""){
              ArrayList<String> old_seq = new ArrayList<String>(Arrays.asList(fasta_map.get(strain).split("")));
              old_seq.remove(index);
              old_seq.add(index,closest);
              StringBuilder sb = new StringBuilder();
              for (String str : old_seq){
                sb.append(str);
              }
              String new_value = sb.toString();
              new_fasta_map.put(strain,new_value);
              g_replacements = g_replacements + 1;
            } else{
              System.out.println("FOUND NO CLOSEST RELATIVE EVER...?");
            }
            index = index + 1;
          } else{
            fasta_map = remove_column(fasta_map,index);
            orig_fasta_map = remove_column(orig_fasta_map,index);
            c_removals = c_removals + 1;
          }
        } 
        else if(!(Arrays.asList(nuc_list).contains(nuc))){
          strain_node = find_node_with_taxon_label(strain,tree);
          HashMap<Character,List<List<List<String>>>> close_rel = 
              weighted_seq_choose(orig_fasta_map,strain_node,rels,index,0);
          Iterator it = close_rel.entrySet().iterator();
          while (it.hasNext()){
            Map.Entry pairs = (Map.Entry)it.next();
            closest = (String)pairs.getKey();
            rels = (List<List<List<String>>>)pairs.getValue();
          }
          if(closest != null){
            ArrayList<String> old_seq = new ArrayList<String>(Arrays.asList(fasta_map.get(strain).split("")));
            old_seq.remove(index);
            old_seq.add(index,closest);
            StringBuilder sb = new StringBuilder();
            for (String str : old_seq){
              sb.append(str);
            }
            String new_value = sb.toString();
            new_fasta_map.put(strain,new_value);
            base_replacements = base_replacements + 1;
          }
          index = index + 1;
        }
        else {
          index = index + 1;
        }
        s_num = s_num + 1;
      }
    }
    return new_fasta_map;
  }
  
  private static void update_fasta(Fasta fasta, HashMap<String,String> fasta_map){
    Iterator it = fasta_map.entrySet().iterator();
    while (it.hasNext()){
      Map.Entry pairs = (Map.Entry)it.next();
      String strain = (String)pairs.getKey();
      String seq = (String)pairs.getValue();
      fasta.setSequence(strain,seq);
    }
  }
}
