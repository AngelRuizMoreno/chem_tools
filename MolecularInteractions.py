from plip.structure.preparation import PDBComplex
import pandas as pd

def run_plip(pdb_file:str=""):
    plip= PDBComplex()
    plip.load_pdb(pdb_file)
    ligands=plip.ligands
    bsids=[f"{lig.hetid}:{lig.chain}:{lig.position}" for lig in ligands]
    plip.analyze()

    def _get_interactions_info(ligand_id:str=""):
    
        interactions = plip.interaction_sets[ligand_id] 
        table=pd.DataFrame()

        for index,interaction in enumerate(interactions.all_itypes):
            interaction_type=str(type(interaction)).split('.')[-1].replace("'>","")
            table.loc[index,'Residue']=interaction.restype+str(interaction.resnr)
            table.loc[index,'Chain']=interaction.reschain
            table.loc[index,'Ligand']=interaction.restype_l+str(interaction.resnr_l)
            
            if interaction_type == 'hbond':
                table.loc[index,'Type']='H-bond'
                table.loc[index,'Acceptor']=interaction.atype
                table.loc[index,'AcceptorIdx']=interaction.a.idx
                table.loc[index,'Donor']=interaction.dtype
                table.loc[index,'DonorIdx']=interaction.d.idx
                table.loc[index,'DistanceAD']=interaction.distance_ad
                table.loc[index,'DistanceAH']=interaction.distance_ah
                table.loc[index,'Angle']=interaction.angle
                table.loc[index,'Force']=interaction.type
                table.loc[index,'ProtIsDon']=interaction.protisdon
            
            elif interaction_type == 'pication':
                table.loc[index,'Type']='Pi-cation'
                table.loc[index,'Charge']=interaction.charge.type
                table.loc[index,'ChargedAtoms']=",".join([i.type for i in interaction.charge.atoms])
                table.loc[index,'Force']=interaction.type
                table.loc[index,'RingType']=interaction.ring.type
                table.loc[index,'RingAtoms']=",".join([i.type for i in interaction.ring.atoms])
                table.loc[index,'RingAtomsIdx']=",".join([str(i.idx) for i in interaction.ring.atoms])
        
            elif interaction_type == 'pistack':
                table.loc[index,'Type']='Pi-stacking'
                table.loc[index,'StackingType']=interaction.type
                table.loc[index,'RecRingType']=interaction.proteinring.type
                table.loc[index,'LigRingType']=interaction.ligandring.type
                table.loc[index,'RecRingAtoms']=",".join([i.type for i in interaction.proteinring.atoms])
                table.loc[index,'RecAtomsIdx']=",".join([str(i.idx) for i in interaction.proteinring.atoms])
                table.loc[index,'LigRingAtoms']=",".join([i.type for i in interaction.ligandring.atoms])
                table.loc[index,'LigRingAtomsIdx']=",".join([str(i.idx) for i in interaction.ligandring.atoms])
                table.loc[index,'Distance']=interaction.distance
                table.loc[index,'Angle']=interaction.angle
                table.loc[index,'Offset']=interaction.offset        
            
            elif interaction_type=='saltbridge':
                table.loc[index,'Type']='Salt-bridge'
                table.loc[index,'NegAtoms']=",".join([i.type for i in interaction.negative.atoms])
                table.loc[index,'NegAtomsIdx']=",".join([str(i.idx) for i in interaction.negative.atoms])
                table.loc[index,'PosAtoms']=",".join([i.type for i in interaction.positive.atoms])
                table.loc[index,'PosAtomsIdx']=",".join([str(i.idx) for i in interaction.positive.atoms])
                table.loc[index,'Distance']=interaction.distance
                table.loc[index,'ProtIsPos']=interaction.protispos
                
            elif interaction_type == 'hydroph_interaction':
                table.loc[index,'Type']='Hydrophobic'
                table.loc[index,'RecAtom']=interaction.bsatom.type
                table.loc[index,'RecAtomIdx']=interaction.bsatom.idx
                table.loc[index,'LigAtom']=interaction.ligatom.type
                table.loc[index,'LigAtomIdx']=interaction.ligatom.idx
                table.loc[index,'Distance']=interaction.distance
                
            elif interaction_type == 'waterbridge':
                table.loc[index,'Type']='Water-bridge'
                table.loc[index,'AccType']=interaction.atype
                table.loc[index,'DonType']=interaction.dtype
                table.loc[index,'WaterIdx']=interaction.water_orig_idx
                table.loc[index,'DistanceAWat']=interaction.distance_aw
                table.loc[index,'DistanceDWat']=interaction.distance_dw
                table.loc[index,'AngleDon']=interaction.d_angle
                table.loc[index,'AngleWat']=interaction.w_angle
                table.loc[index,'ProtIsDon']=interaction.protisdon
        
            elif interaction_type == 'halogenbond':
                table.loc[index,'Type']='X-bond'
                table.loc[index,'Acceptor']=interaction.acctype
                table.loc[index,'Donor']=interaction.donortype
                table.loc[index,'Distance']=interaction.distance
                table.loc[index,'DonAngle']=interaction.don_angle
                table.loc[index,'AccAngle']=interaction.acc_angle
        
            elif interaction_type=='metal_complex':
                table.loc[index,'Type']='Metal-complex'
                table.loc[index,'MetalType']=interaction.metal.type
                table.loc[index,'Idx']=interaction.metal.idx
                table.loc[index,'TargetType']=interaction.target_type
                table.loc[index,'FunctGroup']=interaction.target.fgroup
                table.loc[index,'Geometry']=interaction.geometry
                table.loc[index,'Distance']=interaction.distance
                table.loc[index,'Location']=interaction.location
        
        return table

    for entry in bsids:
        result=_get_interactions_info(entry)
        result.to_csv(f"{pdb_file.replace('.pdb',f'_{entry}.tsv')}",sep="\t",index=False)