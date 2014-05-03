package ecosim;

import java.util.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.IOException;

public class SequenceFix{ 

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
        HashMap<String,String> logSettings = new HashMap<String,String>();
        logSettings.put("Error","Could not find all required settings in 'seq_clean_settings.txt'");
        logSettings.put("-- NEED", "newick file, fasta file, cutoff, new fasta file " +
                        "name, new log name. Each should be entered on a new line " +
                        "in the specified order.");
        logging(logSettings,"SEQ_CLEAN_ERROR_LOG.txt");
        return;
      }} else{
        HashMap<String,String> logFile = new HashMap<String,String>();
        logFile.put("Error","Could not find 'seq_clean_settings.txt'");
        logging(logFile, "SEQ_CLEAN_ERROR_LOG.txt");
        return;
      }
      // load in fasta file and get strain:sequences in a hashmap
      Fasta fasta = new Fasta(fasta_file);
      HashMap<String,String> fasta_map = fasta.getSequenceHash();
      @SuppressWarnings("unchecked")
      HashMap<String,String> fasta_map_copy = (HashMap<String,String>)fasta_map.clone();
      
      HashMap<HashMap<String,String>,HashMap<String,String>> map_log = 
        replacer(newick_file,fasta_map_copy,cutoff);
      Iterator it = map_log.entrySet().iterator();
      while (it.hasNext()){
        Map.Entry pairs = (Map.Entry)it.next();
        @SuppressWarnings("unchecked")
        HashMap<String,String> new_fasta_map = (HashMap<String,String>)pairs.getKey();
        @SuppressWarnings("unchecked")
        HashMap<String,String> logs = (HashMap<String,String>)pairs.getValue();
        System.out.println(new_fasta_map);
        update_fasta(fasta,new_fasta_map);
        fasta.save(fasta_out);
        
        logs.put("Number of strains in fasta file",Integer.toString(new_fasta_map.size()));
        logging(logs,log_name);
      }
      //System.out.println(new_fasta_map);
      //update_fasta(fasta,new_fasta_map);
      //fasta.save(fasta_out);
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
  
  /**
   * Returns the descendants of the given NewickTreeNode. Includes non-leaf nodes.
   * 
   * @param node NewickTreeNode whose descendants are wanted.
   * @return List<NewickTreeNode> containing descendants of the given node.
   */
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
        // Sorting
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
  
  /** Returns the next closest relatives of the given strain
   *  for cases where there are no acceptable closest relatives
   *  returned by closest_relative, i.e. the sister node and its
   *  children. This function goes up to the parent of the
   *  parent of the given strain and then searches down for
   *  relatives.
   *  
   *  @param strain_node NewickTreeNode of data wanted.
   *  @return HashMap<NewickTreeNode,List<List<String>>> containing 
   *          the parent's parent of the node and the 
   *          List of [relatives,distance].
   */
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
    
    // logging method. may or may not need when added to ecosim
    private static void logging(HashMap<String,String> logs, String log_name){
      Date date = new Date();
      String string = "Date: " + date.toString() + "\n" + "\n";
      Iterator it = logs.entrySet().iterator();
      while (it.hasNext()){
        Map.Entry pairs = (Map.Entry)it.next();
        String k = (String)pairs.getKey();
        String logk = (String)pairs.getValue();
        string = string + k + ": " + logk + "\n" + "\n";
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
      
  /** Calculates the percentage of strains that also have a gap
   *  at the same index.
   *  
   *  @param fasta_map HashMap<String,String> containing the strain:sequence.
   *  @param index int index being looked at.
   *  @return percentage of strains that have a gap at the same index. 
   */
  private static double percent_gap(HashMap<String,String> fasta_map, int index){
    int num_strains = fasta_map.size();
    int gaps = 0;
    Iterator it = fasta_map.entrySet().iterator();
    while (it.hasNext()){
      Map.Entry pairs = (Map.Entry)it.next();
      String value = (String)pairs.getValue();
      try{
      if (Character.toString(value.charAt(index)).equals("-")){
          gaps = gaps + 1;
        }} catch(IndexOutOfBoundsException e) {
          gaps = gaps + 1;
        }
    }
    return (double) gaps/num_strains;
  }
  
  /** Removes the nucleotide at the given index for all sequences.
   *  
   *  @param fasta_map HashMap<String,String> containing the strain:sequence.
   *  @param index int index being looked at.
   *  @return HashMap<String,String> containing the strain:new sequence with removed 
   *          nucleotide.
   */
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
  
  /** Chooses the sequence that is most related to the given strain.
   *  
   *  @param strain_fasta_map HashMap<String,String> containing the strain:sequence.
   *  @param strain_node NewickTreeNode containing the strain being looked at.
   *  @param rels List<List<List<String>>> containing the closest relatives of strain_node.
   *  @param index index being looked at.
   *  @param level int representing number of upward requests requested to go further up the tree. 
   *  @return HashMap<String,List<List<List<String>>>> of sequence chosen:other closest relatives
   */
  @SuppressWarnings("unchecked")
  private static HashMap<String,List<List<List<String>>>> weighted_seq_choose(HashMap<String,String> strain_fasta_map, 
                                     NewickTreeNode strain_node, 
                                     List<List<List<String>>> rels, 
                                     int index, int level){
    HashMap<String,Integer> nucs = new HashMap<String,Integer>();
    String[] nuc_list =  new String[] {"A","T","C","G","a","t","c","g"};
    List<List<String>> rel_level = rels.get(level); 
    for(List<String> r : rel_level){
      try{
        String nucleotide = Character.toString((strain_fasta_map.get(r.get(0))).charAt(index));
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
    List<String> pop = new ArrayList<String>();
    Iterator it = nucs.entrySet().iterator();
    while (it.hasNext()){
      Map.Entry pairs = (Map.Entry)it.next();
      String key = (String)pairs.getKey();
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
      String choice = pop.get(rand_num);
      HashMap<String,List<List<List<String>>>> char_rels = new HashMap<String,List<List<List<String>>>>();
      char_rels.put(choice,rels);
      return char_rels;
  }
  
  /** Returns a HashMap<String,String> of strain name:sequence with the sequences
   *  cleaned up (gaps removed/replaced by comparing with the closest relative).
   *  
   *  @param newick_file Name of the newick file (corrosponding to the fasta file).
   *  @param fasta_map HashMap<String,String> containing the strain name:sequence
   *  @return HashMap<String,String> containing the strain name:new sequence
   */
  // CHANGE TO HASHMAP<HASHMAP<STRING,STRING>,HASHMAP<STRING,STRING>> TO RETURN
  // fasta_map and logs
  @SuppressWarnings("unchecked")
  private static HashMap<HashMap<String,String>,HashMap<String,String>> replacer(String newick_file, 
                                                 HashMap<String,String> fasta_map, 
                                                 double cutoff) throws IOException, InvalidNewickException{
    //create tree
    String stree = new Scanner(new File(newick_file)).useDelimiter("\\z").next();
    NewickTree tree = new NewickTree(stree);
    HashMap<String,String> orig_fasta_map = (HashMap<String,String>)fasta_map.clone();
    HashMap<String,String> log = new HashMap<String,String>();
    String[] nuc_list =  new String[] {"A","T","C","G","a","t","c","g"};
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
      NewickTreeNode strain_node = null;
      for(int i = 0; i < sequence_copy.length(); i++){
        String nuc = Character.toString(sequence_copy.charAt(i));
        if(nuc.equals("-")){
          double p_gaps = percent_gap(fasta_map,index);
          if(p_gaps <= cutoff){
            if(strain_node == null){
              strain_node = new NewickTreeNode();
              strain_node = find_node_with_taxon_label(strain,tree);
            }
            
            HashMap<String,List<List<List<String>>>> close_rel = 
              weighted_seq_choose(orig_fasta_map,strain_node,rels,index,0);
            Iterator it = close_rel.entrySet().iterator();
            while (it.hasNext()){
              Map.Entry pairs = (Map.Entry)it.next();
              closest = (String)pairs.getKey();
              rels = (List<List<List<String>>>)pairs.getValue();
            }
            if(!(closest.equals(""))){
              ArrayList<String> old_seq = new ArrayList<String>(Arrays.asList(fasta_map.get(strain).split("")));
              old_seq.remove(0);
              old_seq.remove(index);
              old_seq.add(index,closest);
              StringBuilder sb = new StringBuilder();
              for (String str : old_seq){
                sb.append(str);
              }
              String new_value = sb.toString();
              fasta_map.remove(strain);
              fasta_map.put(strain,new_value);
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
        else if((Arrays.asList(nuc_list).contains(nuc)) == false){
          strain_node = find_node_with_taxon_label(strain,tree);
          HashMap<String,List<List<List<String>>>> close_rel = 
              weighted_seq_choose(orig_fasta_map,strain_node,rels,index,0);
          Iterator it0 = close_rel.entrySet().iterator();
          while (it0.hasNext()){
            Map.Entry pairs = (Map.Entry)it0.next();
            closest = (String)pairs.getKey();
            rels = (List<List<List<String>>>)pairs.getValue();
          }
          if(!(closest.equals(""))){
            ArrayList<String> old_seq = new ArrayList<String>(Arrays.asList(fasta_map.get(strain).split("")));
            // always gives back empty first value, so remove it.
            old_seq.remove(0);
            old_seq.remove(index);
            old_seq.add(index,closest);
            StringBuilder sb = new StringBuilder();
            for (String str : old_seq){
              sb.append(str);
            }
            String new_value = sb.toString();
            fasta_map.remove(strain);
            fasta_map.put(strain,new_value);
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
    HashMap<HashMap<String,String>,HashMap<String,String>> map_log = 
      new HashMap<HashMap<String,String>,HashMap<String,String>>();
    log.put("Strins analyzed",Arrays.toString(fasta_map.keySet().toArray()));
    log.put("Cutoff percentage",Double.toString(cutoff));
    log.put("Number of strains in Newick Tree", Integer.toString(fasta_map.keySet().size()));
    log.put("Number of column removals", Integer.toString(c_removals));
    log.put("Number of gap replacements", Integer.toString(g_replacements));
    log.put("Number of irregular base replacements", Integer.toString(base_replacements));
    map_log.put(fasta_map,log);
    return map_log;
  }
  
  /** Updates the given fasta object with the new sequences.
   *  
   *  @param fasta Fasta object 
   *  @param fasta_map HashMap<String,String> containing the strain:new sequence
   */
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
