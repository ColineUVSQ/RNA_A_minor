'''
Created on 10 avr. 2019

@author: coline
'''

EXTENSION_PATH = 'Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles_new_dernier/result_graphes_extension_%s/'
EXTENSION_PATH_TAILLE = 'Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles_new_dernier/graphes_extension_taille_%s/' 
EXTENSION_TOUTES_ARETES = 'Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles_new_dernier/'

NEW_EXTENSION_PATH_TAILLE = 'Nouvelles_donnees/'

PATH_MMCIF = "Fichiers_mmcif/"

GROUPES_TOUTES_ARETES_MAX_4_10_07 = [['4V88_A5_25_47', '5DM6_X_25_34', '5FDU_1A_25_78', '5J7L_DA_25_10', '3JCS_1_25_46', '4YAZ_R_36_25'],
                                     ['1FJG_A_48_17', '5J5B_BA_48_14', '5DM6_X_25_15'],
                                     ['1U9S_A_58_11', '5FDU_1A_48_19', '5J7L_DA_30_15', '5FDU_1A_30_17', '5J5B_BA_48_23', '1FJG_A_48_8', '5J7L_DA_48_15', '4V9F_0_30_4'],
                                     ['5DM6_X_328_2', '5J7L_DA_62_5', '4V9F_0_62_12', '5FDU_1A_62_14', '3JCS_1_25_16'],
                                     ['4V9F_0_48_13', '5J7L_DA_48_20', '5J7L_DA_48_1', '5FDU_1A_134_3', '5J7L_DA_134_1', '5DM6_X_134_2', '4V9F_0_127_6', '4V88_A5_48_3', '5J7L_DA_48_30', '5DM6_X_48_10', '1FJG_A_271_1', '5FDU_1A_137_6', '4V9F_0_137_5', '5DM6_X_48_28', '3JCS_1_48_18', '4V9F_0_48_16', '5J5B_BA_48_5', '5FDU_1A_48_25']]

GROUPE_GNRA = ['1U9S_A_58_11', '3UCZ_R_62_15', '4V9F_0_44_3', '1FJG_A_48_11', '1FJG_A_58_23', '4V9F_0_48_26', '5J7L_DA_50_21']
GROUPE_GNRA_ETENDU = ['4V88_A5_25_47', '5DM6_X_25_34', '5FDU_1A_25_78', '5J7L_DA_25_10', '3JCS_1_25_46', '4YAZ_R_36_25', '5FDU_1A_48_19', '5J7L_DA_30_15', '5FDU_1A_30_17', '5J5B_BA_48_23', '1FJG_A_48_8', '5J7L_DA_48_15', '4V9F_0_30_4', '5FDU_1A_25_68', '4V88_A5_25_30', '4V9F_0_25_56', '5J7L_DA_25_12'] ## boucle de 4 nts de type GNRA mais pas la liaison tertiaire non canonique
GROUPE_GBULGE = ['4FAU_A_207_4', '5FDU_1A_272_1', '5J7L_DA_272_2', '4V9F_0_207_3']
GROUPE_ARICH = ['5FDU_1A_137_6', '5DM6_X_48_10', '5FDU_1A_134_3', '5DM6_X_134_2', '4V88_A5_48_3', '5J7L_DA_134_1', '4V9F_0_137_5', '4V9F_0_134_5', '5J7L_DA_48_1', '5FDU_1A_74_7', '4V9F_0_127_6']
GROUPE_ARICH_ETENDU = ['3JCS_1_48_18', '1FJG_A_271_1', '4V9F_0_48_2', '4V9F_0_48_16', '5J7L_DA_48_30'] # manque la liaison THW

GROUPE_ARICH_DE_07 = ['5FDU_1A_137_6', '5DM6_X_48_10', '5FDU_1A_134_3', '5DM6_X_134_2', '4V88_A5_48_3', '5J7L_DA_134_1', '4V9F_0_137_5', '5J7L_DA_48_1', '4V9F_0_127_6']
GROUPE_DENSE_DE_07 = ['5FDU_1A_137_6', '5DM6_X_48_10', '5FDU_1A_134_3', '5DM6_X_134_2', '4V88_A5_48_3', '5J7L_DA_134_1', '4V9F_0_137_5', '5J7L_DA_48_1', '4V9F_0_127_6', '5J7L_DA_48_20', '5FDU_1A_48_25', '5DM6_X_48_28', '4V9F_0_48_13', '1FJG_A_271_1']


GROUPE_ETOILE_GNRA = ['1FJG_A_58_23', '5J5B_BA_58_3', '4V9F_0_44_3', '4V88_A6_17_55']

HOMOLOGUES = [['4V9F_0_62_12', '5J7L_DA_62_5', '5DM6_X_328_2', '5FDU_1A_62_14', '3JCS_1_25_16'],
                  ['4V9F_0_30_4', '5J7L_DA_30_15', '5FDU_1A_30_17'],
                  ['4V9F_0_48_13', '5J7L_DA_48_20', '5DM6_X_48_28', '5FDU_1A_48_25'],
                  ['4V9F_0_207_3', '5J7L_DA_272_2', '5FDU_1A_272_1'],
                  ['4V9F_0_25_56', '5J7L_DA_25_12', '5FDU_1A_25_68', '4V88_A5_25_30'],
                  ['4V9F_0_137_5', '5J7L_DA_48_1', '5FDU_1A_137_6', '5DM6_X_48_10', '4V88_A5_48_3', '3JCS_1_48_18'],
                  ['4V9F_0_134_5', '5FDU_1A_74_7', '5DM6_X_127_7'],
                  ['4V9F_0_48_21', '5J7L_DA_197_4', '5FDU_1A_197_3', '5DM6_X_48_9'],
                  ['4V9F_0_127_6', '5J7L_DA_134_1', '5FDU_1A_134_3', '5DM6_X_134_2'],
                  ['4V9F_0_287_2', '5DM6_X_25_15'],
                  ['5J7L_DA_25_10', '5DM6_X_25_34', '4V88_A5_25_47', '5FDU_1A_25_78', '3JCS_1_25_46'],
                  ['1FJG_A_48_17', '5J5B_BA_48_14'],
                  ['1FJG_A_48_8', '5J5B_BA_48_23'],
                  ['1FJG_A_294_1', '5J5B_BA_294_2'],
                  ['1FJG_A_58_23', '5J5B_BA_58_3'],
                  ['1FJG_A_138_3', '5J5B_BA_138_2'],
                  ['5J5B_BA_48_7', '4V88_A6_48_12'],
                  ['4V9F_0_44_3', '5J7L_DA_50_21']]


CLUSTERING_PEREZ_VERSION_NON_CAN_2 = [['5J5B_BA_48_14', '1FJG_A_48_17'], 
                                      ['5J7L_DA_25_10', '4V88_A5_25_47', '5FDU_1A_25_78', '5DM6_X_25_34', '3JCS_1_25_46'], 
                                      ['5FDU_1A_134_3', '5J7L_DA_134_1', '5DM6_X_134_2', '4V9F_0_127_6'], 
                                      ['4V88_A6_48_12', '5J5B_BA_48_7'], 
                                      ['4V9F_0_48_13', '5DM6_X_48_28', '5FDU_1A_48_25', '5J7L_DA_48_20'], 
                                      ['5DM6_X_227_2', '5J7L_DA_48_30'], 
                                      ['5FDU_1A_25_68', '4V88_A5_25_30', '5J7L_DA_25_12', '4V9F_0_25_56'], 
                                      ['5J5B_BA_138_2', '1FJG_A_138_3'], 
                                      ['1FJG_A_294_1', '5J5B_BA_294_2'], 
                                      ['4V9F_0_62_12', '5DM6_X_328_2', '5J7L_DA_62_5', '5FDU_1A_62_14', '3JCS_1_25_16'], 
                                      ['5J7L_DA_197_4', '5DM6_X_48_9', '4V9F_0_48_21', '5FDU_1A_197_3'], 
                                      ['5J7L_DA_30_15', '4V9F_0_30_4', '5J5B_BA_48_23', '5FDU_1A_30_17', '1FJG_A_48_8'], 
                                      ['5J7L_DA_48_15', '1FJG_A_48_8', '1U9S_A_58_11', '5FDU_1A_48_19', '5J5B_BA_48_23'], 
                                      ['1FJG_A_271_1', '5DM6_X_48_10'], 
                                      ['5FDU_1A_137_6', '4V88_A5_48_3', '5DM6_X_48_10', '5J7L_DA_48_1', '4V9F_0_137_5', '3JCS_1_48_18'], 
                                      ['5J7L_DA_272_2', '5FDU_1A_272_1'], 
                                      ['4V88_A6_17_55', '1FJG_A_58_23'], 
                                      ['3UCZ_R_62_15', '1U9S_A_58_11'], 
                                      ['4V9F_0_207_3', '5FDU_1A_272_1'], 
                                      ['4V9F_0_44_3', '1FJG_A_58_23'], 
                                      ['5DM6_X_127_7', '5FDU_1A_74_7'], 
                                      ['4V9F_0_48_16', '4V88_A5_48_3'], 
                                      ['5J5B_BA_58_3', '1FJG_A_58_23'], 
                                      ['4V9F_0_134_5', '5FDU_1A_74_7']]

GROUPE_JAUNE_VERT_BLEU = ['5J7L_DA_30_15', '4V9F_0_30_4', '5J5B_BA_48_23', '5FDU_1A_30_17', '1FJG_A_48_8', '5J7L_DA_48_15', '1U9S_A_58_11', '5FDU_1A_48_19']


''' 11/09/19 liste des problemes de doubles liaisons '''
liste_pbs = ['1ond_2', '1ond_5', '4v7s_2', '4v7s_6', '4v7s_20', '4v7s_57', '4wqy_6', '4wqy_14', '4wqy_23', '4wqy_54', '5mdv_3', '5j7l_20', '5j7l_51', '3fwo_9', '5wit_40', '5wit_50', '5wit_51', '1hnx_1', '1hnx_9', '5t7v_1', '5t7v_8', '5t7v_12', '4v5c_7', '4v5c_21', '4v5c_45', '4v5c_47', '3j9w_2', '1sm1_3', '1sm1_4', '1sm1_8', '4u52_8', '4ox9_1', '1hr2_1', '1hr2_2', '5e7k_6', '5e7k_13', '5e7k_28', '5e7k_44', '5e7k_46', '5e7k_47', '6b4v_13', '6b4v_32', '6b4v_43', '1vq6_3', '1vq6_19', '4ujc_2', '4ujc_3', '4dv4_5', '5jvg_2', '5jvg_11', '4io9_2', '4io9_12', '4io9_15', '5ib8_6', '5ib8_31', '5ib8_39', '5ib8_49', '5ib8_51', '5ib8_52', '4v4a_4', '4v4a_7', '4v4a_9', '1y69_12', '1y69_13', '5mdz_3', '5mdz_37', '6bz8_12', '6bz8_15', '6bz8_23', '6bz8_28', '6bz8_41', '6bz8_55', '4wqu_3', '4wqu_7', '4wqu_13', '4wqu_22', '4wqu_55', '6n8o_7', '5dfe_4', '5dfe_10', '5dfe_11', '5dfe_14', '4v6i_3', '1yjn_17', '1yjn_20', '5kcs_2', '5kcs_3', '5kcs_17', '4v4y_2', '6q8y_10', '4dr3_1', '4dr3_8', '6fkr_14', '6fkr_24', '6fkr_28', '6fkr_43', '6fkr_51', '6cfj_12', '6cfj_39', '6cfj_48', '6cfj_49', '1njo_1', '1njo_2', '1njo_4', '1njo_15', '5mgp_11', '5mgp_15', '4v7k_1', '4v7k_7', '4v7k_21', '5njt_5', '5njt_18', '1xmo_11', '4wpo_3', '4wpo_8', '4wpo_18', '4wpo_37', '5f8k_31', '5f8k_34', '5f8k_36', '5f8k_50', '3cc4_15', '3cc4_21', '4v6e_8', '4v6e_18', '4v6e_51', '6c5l_1', '6c5l_10', '6c5l_32', '6c5l_44', '6c5l_46', '6c5l_47', '6czr_5', '6czr_8', '6czr_37', '5gag_12', '6gz3_7', '3j6x_5', '4woi_11', '4woi_16', '4woi_44', '4woi_56', '4woi_61', '5lzf_8', '1nwy_1', '1nwy_4', '4v9r_26', '4v9r_31', '4v9r_34', '4v9r_43', '4v9r_53', '5a2q_3', '4v67_24', '4v67_33', '4v67_46', '4v67_58', '5gak_1', '5hcp_35', '5hcp_38', '5hcp_41', '5hcp_44', '5hcp_53', '5lmv_2', '2zjr_2', '2zjr_9', '2zjr_12', '4v79_1', '5fdv_38', '5fdv_42', '5fdv_44', '5fdv_57', '5fdv_60', '1vy4_19', '1vy4_38', '1vy4_40', '1vy4_59', '1vy4_63', '3j7z_4', '3j7z_14', '1vqp_21', '4u4z_8', '3ow2_2', '3ow2_10', '5t62_6', '4v8p_10', '4v8p_12', '4v8p_27', '4v8p_30', '4v8p_38', '4v51_9', '4v51_37', '4v51_39', '4v51_50', '4v51_52', '4v51_58', '4v9f_1', '4k0k_4', '4k0k_12', '2otl_4', '2otl_21', '4v8h_1', '4v8h_26', '4v8h_29', '4v8h_37', '2uu9_13', '4tue_1', '4tue_9', '4tue_10', '4tue_14', '6h4n_14', '4v8d_16', '4v8d_29', '4v8d_42', '4v8d_51', '4u4n_1', '4u4n_20', '5w4k_39', '5w4k_52', '5w4k_53', '4v9j_1', '4v9j_39', '3ccr_19', '3jce_12', '3jce_13', '1jzy_5', '1jzy_9', '1jzy_12', '4yzv_2', '4yzv_8', '4yzv_23', '4yzv_34', '4yzv_53', '6d8o_1', '6d8o_2', '5nrg_2', '2ogo_1', '2ogo_4', '2ogo_6', '5jc9_22', '5jc9_49', '1vy6_24', '1vy6_30', '1vy6_33', '1vy6_47', '4u3u_1', '4u3u_6', '4u3u_12', '4u3u_22', '4wq1_8', '4wq1_20', '4wq1_30', '4wq1_36', '3i55_21', '5uyp_11', '2zjp_2', '2zjp_9', '6hcm_5', '5hcr_32', '5hcr_33', '5hcr_36', '5hcr_39', '5hcr_51', '1k73_23', '5gae_1', '5gae_15', '5uq8_22', '4v69_6', '5hl7_2', '4b3t_10', '5lzd_16', '5j8a_21', '5j8a_52', '4lfz_2', '4lfz_22', '4lfz_26', '4lfz_35', '1kc8_21', '3pip_1', '3pip_7', '4wt8_5', '4wt8_38', '4wt8_42', '4wt8_43', '4v9h_19', '5lyb_16', '6d8m_1', '6d8m_2', '2ogm_3', '2ogm_6', '2ogm_12', '4wf9_2', '4v8f_17', '4v8f_21', '4v8f_24', '4v8f_39', '4v8f_40', '4v8f_45', '3iin_1', '6g5i_2', '5ndw_1', '4duz_6', '4duz_8', '1ffk_17', '1ffk_19', '4v8j_1', '4v8j_35', '4v8j_50', '4v8j_52', '4v8j_53', '4u3m_1', '5j01_1', '4v9d_5', '4v9d_17', '4v9d_19', '4v9d_32', '4v9d_33', '4v9d_46', '4v9d_55', '4v53_17', '4v53_29', '4v53_38', '4v53_44', '2uxb_11', '6boh_7', '6boh_49', '6boh_52', '3cme_15', '4lnt_2', '4lnt_25', '4lnt_36', '4lnt_61', '4adv_4', '1vq4_3', '1vq4_19', '4v84_1', '4v84_44', '3jq4_7', '4u50_1', '4u50_8', '6hhq_1', '6hhq_7', '6hhq_9', '1nkw_2', '1nkw_5', '6q9a_2', '1yij_15', '3j9y_14', '4v5m_11', '4v5m_20', '1hnz_13', '4v5a_1', '4v5a_7', '4v5a_26', '4v5a_31', '4v5a_33', '4v5a_38', '5mrf_8', '5mrf_14', '5xxb_8', '5xxb_10', '1vq8_18', '2o45_1', '4e8t_2', '6nd6_1', '6nd6_13', '6nd6_38', '6nd6_49', '6nd6_50', '1j5e_1', '4v88_23', '6bz6_11', '6bz6_14', '4v5y_37', '4v5y_45', '4xej_15', '4xej_25', '4xej_28', '4xej_40', '4v6g_33', '4v6g_36', '4v6g_44', '5j91_23', '5j91_53', '6cae_13', '6cae_43', '6cae_54', '6cae_55', '6ddd_13', '4v4w_4', '1gid_1', '1gid_2', '4dr1_7', '4dr1_8', '6mtb_2', '1njm_2', '1njm_4', '1njm_15', '4v7i_1', '4v7i_7', '5dge_8', '6gc8_8', '1n33_1', '1n33_8', '5j4d_7', '5j4d_15', '5j4d_20', '5j4d_50', '4v7e_7', '6n8m_5', '4y4p_44', '4y4p_56', '4y4p_57', '4v6k_10', '4v6k_13', '4v6k_14', '6gsk_6', '6gsk_23', '6gsk_37', '6gsk_40', '4tuc_10', '4tuc_11', '4tuc_15', '1vqn_19', '4v49_3', '4v49_6', '4v49_8', '1l8v_1', '1l8v_2', '4v8n_10', '4v8n_11', '4v8n_22', '4v8n_38', '4ioa_2', '4ioa_10', '4ioa_13', '4wf1_6', '4wf1_8', '4wf1_17', '4wf1_50', '2otj_14', '2otj_20', '5uyl_11', '5uyl_22', '5myj_11', '1q81_3', '1q81_18', '1q81_19', '6awd_2', '6awd_4', '3cpw_4', '3cpw_20', '5h5u_6', '5h5u_22', '4v57_14', '4v57_38', '4v57_45', '4v57_51', '4v57_55', '2r8s_1', '6hcq_2', '5lzt_1', '5lzt_8', '5lzt_11', '1zzn_1', '4u1u_5', '4u1u_8', '4u1u_50', '4nxm_4', '4nxm_6', '6dzp_14', '5lzx_13', '6q97_4', '4v9l_4', '6eml_5', '6i7v_7', '6i7v_13', '6i7v_33', '6i7v_42', '4uy8_8', '4uy8_19', '6gwt_9', '5ngm_13', '4v8b_16', '4v8b_28', '4v8b_41', '4v8b_42', '4v8b_50', '6ftj_7', '6ftj_9', '5fci_4', '3ccl_14', '6gz5_5', '6gz5_7', '6gz5_11', '6gz5_13', '3jbu_11', '3jbu_12', '3jbu_16', '3jbu_19', '5vpo_8', '5vpo_9', '5vpo_11', '1kd1_3', '1kd1_21', '4v73_8', '4v73_11', '2vqe_9', '5ndk_25', '5ndk_33', '5ndk_40', '5ndk_55', '3j7p_1', '3j7p_3', '3j7p_18', '4csu_1', '4lel_2', '4lel_10', '4lel_30', '4wzo_6', '4wzo_27', '4wzo_38', '4wzo_41', '4wzo_42', '4wzo_44', '5ndg_6', '5ndg_15', '5lks_1', '5lks_2', '5lks_8', '3jag_1', '3jag_14', '5apn_5', '3j5l_4', '3j5l_13', '1qvg_3', '1qvg_18', '4v61_8', '5d8b_6', '5d8b_13', '5d8b_17', '5d8b_21', '6enf_15', '5nd9_2', '5o61_11', '4v7a_10', '5el7_7', '5el7_16', '5el7_31', '5el7_47', '5el7_50', '5el7_52', '5el7_55', '4v6o_1', '4v5q_1', '4v5q_29', '4v5q_45', '6gb2_1', '6fyy_2', '6fyy_3', '6fyy_4', '5kpv_23', '5j5b_23', '5j5b_53', '5jte_7', '4v6c_6', '4v6c_56', '3cc2_16', '5it8_19', '5it8_48', '6cfl_35', '6cfl_38', '6cfl_41', '6cfl_52', '6cfl_55', '1nji_22', '4v7m_24', '4v7m_30', '4v7m_34', '4v7m_51', '4v4s_4', '4v4s_5', '4v4s_6', '4ypb_32', '4ypb_40', '4ypb_50', '6gc0_5', '5vyc_4', '5vyc_29', '4v5e_10', '4v5e_25', '4v5e_34', '4v5e_38', '4v5e_42', '4v5e_45', '6h58_7', '6h58_30', '4u25_4', '4u25_6', '4u25_14', '4u25_28', '4u25_43', '4u25_44', '4u25_50', '4d67_1', '4d67_3', '4ug0_9', '4ug0_17', '4w2h_15', '4w2h_26', '4w2h_32', '4w2h_35', '4w2h_50', '5hd1_36', '5hd1_39', '5hd1_43', '5hd1_55', '6fxc_10', '6fxc_27', '1j5a_5', '1j5a_9', '1j5a_12', '4lf8_8', '5jus_5', '4wra_16', '4wra_21', '4wra_26', '4wra_41', '4wra_46', '4wra_51', '4v7u_6', '4v7u_18', '4v7u_25', '4v7u_48', '4v7u_54', '5t2a_1', '5t6r_6', '4w4g_2', '4w4g_32', '4w4g_50', '4dv2_6', '6gxo_7', '6gxo_12', '4a2i_3', '3cma_17', '3cma_24', '4v7y_33', '4v7y_35', '4v7y_49', '4v7y_56', '6buw_9', '6buw_10', '6buw_12', '4uje_2', '4uje_3', '5j3c_1', '5j3c_11', '5j3c_12', '5j3c_14', '4v4q_16', '4v4q_32', '4v4q_40', '4v4q_47', '4v6a_4', '4v6a_11', '4v6a_23', '4v6a_27', '4v6a_29', '4v6a_46', '1m90_3', '1m90_23', '5dm7_1', '5dm7_3', '2qa4_3', '2qa4_16', '5kpx_6', '5kpx_24', '4v6m_16', '4v5s_34', '4v5s_52', '4v5s_55', '6n8k_1', '6n8k_3', '4v7c_7', '5j4b_38', '5j4b_41', '5j4b_44', '5j4b_57', '6n9e_29', '6n9e_31', '5el5_9', '5el5_27', '5el5_38', '5el5_42', '5hkv_4', '4v5k_1', '4v5k_6', '4v5k_24', '5dat_6', '5aka_3', '4v6u_11', '4v6u_14', '4w2f_37', '4w2f_39', '4w2f_54', '4w2f_57', '4u56_1', '4u56_8', '4dv0_6', '6gxm_11', '5t2c_2', '5t2c_6', '5tcu_2', '4faq_1', '4wro_6', '4wro_20', '4wro_28', '4wro_30', '4wro_33', '4wro_50', '6ft6_3', '4v4i_10', '3cf5_5', '3cf5_11', '3cf5_14', '4e8r_2', '6ha1_22', '4v7w_7', '4v7w_31', '4v7w_33', '4v7w_44', '4v7w_46', '4v7w_53', '2o43_6', '4u27_5', '4u27_7', '4u27_16', '4u27_52', '1jj2_21', '3cxc_18', '4v90_18', '4v90_24', '4v5g_6', '4v5g_15', '4v5g_27', '4v5g_28', '4v5g_29', '4v5g_41', '5tbw_11', '5tbw_36', '6gzq_1', '4wsm_17', '4wsm_25', '4wsm_31', '4wsm_41', '4wsm_43', '5jb3_4', '3jbo_2', '2a2e_1', '6i0y_4', '6i0y_14', '6q95_27', '3ccv_14', '4v9n_2', '4v9n_41', '4v9n_46', '4v9n_60', '5lzz_1', '4d5y_1', '4d5y_3', '5uyn_2', '5uyn_10', '1yi2_15', '4v55_32', '4v55_35', '4v55_41', '5wfk_15', '5lzv_11', '4v9b_4', '4v9b_5', '4v9b_28', '4v9b_34', '4v9b_49', '4tua_9', '4tua_10', '4tua_15', '1vql_20', '5u9g_12', '4ioc_1', '4ioc_8', '4ioc_9', '4ioc_12', '5ibb_7', '5ibb_31', '5ibb_48', '5ibb_50', '5ibb_52', '5wfs_13', '1xbp_4', '1xbp_11', '4v63_17', '4v63_22', '4v63_27', '4v63_29', '4v63_30', '4v63_40', '1z58_3', '1z58_5', '4v8x_1', '4v8x_8', '4v8x_17', '4v8x_31', '4v8x_32', '4v8x_46', '6o97_15', '6o97_40', '6o97_50', '6o97_51', '4u4r_1', '4u4r_3', '4u4r_8', '4u4r_22', '3j7r_1', '3j7r_2', '3j7r_11', '3j7r_17', '3jai_12', '5doy_13', '5doy_38', '5doy_48', '5doy_49', '5o2r_16', '5lzb_8', '4b3r_11', '4tud_1', '4tud_10', '4tud_11', '4tud_15', '4tud_43', '6gsl_9', '6gsl_12', '6gsl_20', '6gsl_32', '6gsl_37', '6gsl_50', '4u3n_8', '4v8i_26', '4v8i_31', '4v8i_33', '4v8i_51', '1q86_14', '5uyk_8', '5uyk_12', '5uyk_17', '6bok_6', '6bok_10', '6bok_13', '6bok_26', '6bok_27', '6bok_53', '5lzs_1', '5lzs_9', '5lzs_14', '4v50_19', '4v50_39', '4v50_44', '5lmo_5', '6enu_18', '2ogn_5', '1jzx_5', '1jzx_9', '1jzx_12', '4p70_1', '4p70_10', '4p70_28', '4p70_31', '4p70_46', '4zer_33', '4zer_35', '4zer_39', '4zer_41', '4zer_44', '4zer_51', '4zer_53', '4v9k_1', '4v9k_11', '4v9k_14', '4v9k_37', '4ji4_6', '1k8a_3', '6ndk_9', '6ndk_10', '6ndk_37', '5hau_35', '5hau_37', '5hau_39', '5hau_43', '5hau_45', '5hau_55', '5hau_56', '3j7o_1', '3j7o_3', '4duy_7', '4u4o_1', '4u4o_23', '5v93_3', '4v8e_4', '5vpp_4', '5vpp_6', '5vpp_7', '5vpp_9', '5vpp_11', '5vpp_18', '5vpp_31', '5vpp_45', '1nwx_1', '1nwx_4', '4v9s_24', '4v9s_28', '4v9s_31', '4v9s_50', '1y0q_1', '5gaf_10', '4wzd_30', '4wzd_42', '4wzd_46', '4wzd_48', '4v74_1', '5j8b_4', '5j8b_30', '2rkj_2', '2rkj_3', '2rkj_7', '2rkj_8', '4v8q_15', '1vy5_35', '1vy5_37', '1vy5_53', '1vy5_58', '5hcq_36', '5hcq_37', '5hcq_41', '5hcq_45', '5hcq_47', '5hcq_57', '1k01_5', '1k01_9', '1k01_12', '2aar_2', '2aar_5', '2aar_11', '2aar_13', '5oql_2', '3j2f_3', '3i56_14', '3i56_19', '5uq7_19', '3dll_2', '3dll_3', '3dll_13', '4khp_4', '4y1n_1', '4y1n_2', '6gxp_11', '3j28_1', '3j28_5', '5kcr_4', '5we6_12', '5we6_14', '6of1_8', '6of1_13', '6of1_42', '6of1_55', '6of1_56', '6n8n_2', '6n8n_6', '6n8n_7', '4v6d_6', '4v6d_21', '4v6d_32', '4v6d_37', '4v6d_56', '5dgf_6', '2ykr_2', '1njn_3', '1njn_14', '4v7j_1', '4v7j_8', '4v7j_21', '4v7j_35', '6cfk_34', '6cfk_35', '6cfk_38', '6cfk_41', '6cfk_52', '6cfk_55', '4dr2_7', '4dr2_10', '5x8p_18', '4v4t_5', '4v4t_6', '4v4t_9', '4l71_2', '4l71_29', '4l71_36', '4l71_38', '4l71_47', '6ddg_14', '5mre_8', '5mre_14', '4v95_14', '4v95_27', '4v95_33', '4v95_45', '4v95_53', '4v5b_6', '4v5b_7', '4v5b_16', '4v5b_23', '4v5b_31', '4v5b_32', '1u6b_1', '5mdw_3', '5mdw_36', '5jbh_5', '6nd5_7', '6nd5_41', '6nd5_50', '6nd5_51', '4l47_2', '4l47_22', '4l47_32', '4l47_53', '3g6e_15', '4v87_13', '4v87_17', '4v87_21', '4v87_46', '4dv5_6', '5ot7_5', '5ot7_8', '1vq7_3', '1vq7_19', '6ha8_20', '5aa0_1', '5wf0_15', '4aqy_10', '4aqy_12', '3j9z_8', '4v5n_12', '4v5n_22', '4wsd_7', '4wsd_20', '4wsd_29', '4wsd_31', '4wsd_34', '4wsd_50', '4www_30', '4www_37', '1yjw_14', '4u53_1', '4u53_11', '6gzx_1', '6gzx_3', '6by1_11', '6by1_15', '6by1_24', '6by1_32', '4v7h_3', '4e8m_1', '1yhq_15', '5xyu_5', '5kps_1', '4v6f_17', '4v6f_47', '4v6f_63', '1xnr_1', '1xnr_12', '3cc7_15', '5a9z_1', '5we4_17', '3cd6_12', '6n8l_1', '4u67_1', '4v7d_3', '4v7d_13', '5anc_7', '4v4z_8', '4v5l_10', '4v5l_15', '4v5l_23', '4v5l_31', '4u51_6', '4z8c_36', '4z8c_37', '4z8c_41', '4z8c_44', '4z8c_56', '6gzz_1', '6gzz_2', '6gzz_3', '6gzz_4', '1hnw_12', '5my1_5', '4v4b_18', '5mdy_3', '4v85_4', '4v85_7', '6d90_7', '6d90_9', '6d90_10', '1vq5_3', '1vq5_19', '5tga_6', '5tga_11', '5tga_14', '5tga_24', '5xym_12', '4v89_3', '4v89_5', '4v89_19', '4v89_26', '5jvh_2', '5jvh_7', '5jvh_12', '5ib7_9', '5ib7_35', '5ib7_41', '5ib7_52', '5ib7_54', '5ib7_56', '4v4n_11', '4v4n_12', '6bz7_8', '6bz7_9', '6bz7_12', '5tgm_14', '5tgm_21', '5tgm_22', '5tgm_26', '5tgm_30', '2o44_1', '2o44_5', '2o44_6', '2o44_7', '2o44_9', '2o44_11', '1vq9_19', '6hma_14', '4v97_31', '4v97_32', '4v97_40', '4v97_49', '4v97_52', '4u20_6', '4u20_9', '4u20_12', '4u20_15', '4u20_31', '4u20_52', '1m1k_21', '6bjx_1', '5ndv_1', '5ndv_6', '5ndv_8', '4lfc_5', '4lfc_6', '5nwy_32', '5xy3_2', '5xy3_7', '4v8g_26', '4v8g_32', '4v8g_34', '4v8g_52', '1vvj_2', '1vvj_10', '1vvj_14', '1vvj_29', '1vvj_37', '1vvj_47', '1k9m_3', '1k9m_20', '6qul_11', '1jzz_5', '1jzz_9', '1jzz_12', '6d8l_1', '6d8l_2', '4v9i_1', '4v9i_32', '4v9i_47', '4v9i_52', '3ccq_15', '5mmm_15', '2uxc_1', '3jcj_17', '3jcj_19', '4v52_39', '4v52_45', '4wr6_4', '4wr6_7', '4wr6_11', '4wr6_18', '4wr6_27', '4wr6_30', '4wr6_33', '4wr6_44', '4wr6_46', '1vqk_21', '4zsn_31', '4zsn_50', '5mlc_12', '3bo2_1', '5me1_5', '2zjq_2', '2zjq_7', '2zjq_10', '5on6_6', '5on6_11', '3cce_12', '5gah_12', '4v64_40', '4v64_45', '6dnc_1', '6dnc_9', '6gaw_3', '5nco_12', '5iqr_3', '5iqr_15', '1vy7_1', '1vy7_13', '1vy7_26', '1vy7_32', '1vy7_35', '1vy7_48', '5fdu_38', '5fdu_42', '5fdu_44', '5fdu_57', '3jbp_2', '3jbp_5', '3jan_6', '3jan_7', '3jan_8', '3jan_9', '3jan_12', '3jan_17', '6bu8_9', '4u4u_1', '4u4u_9', '5obm_1', '5obm_2', '5obm_5', '5obm_8', '4v9q_7', '4v9q_32', '4z3s_8', '4z3s_45', '4z3s_57', '4z3s_58', '3j92_7', '3j92_11', '5lze_12', '3j2h_5', '4v68_1', '4v68_15', '4v68_17', '5gad_12', '1njp_1', '1njp_2', '1njp_4', '1njp_15', '4v7t_2', '4v7t_6', '4v7t_16', '4v7t_18', '4v7t_54', '4p6f_8', '4p6f_9', '4p6f_14', '4p6f_31', '5xyi_1', '4wfb_2', '4w2i_30', '4w2i_38', '4w2i_40', '4w2i_55', '4w2i_59', '4lt8_2', '4lt8_35', '4lt8_55', '5wis_1', '5wis_33', '5wis_36', '5wis_39', '4v5d_9', '4v5d_19', '4v5d_47', '4v5d_49', '5mrc_8', '5mrc_14', '4u24_6', '4u24_8', '4u24_12', '4u24_52', '4u55_9', '4w2e_29', '4v6v_4', '6qkl_3', '6qkl_6', '5v8i_9', '5v8i_20', '5v8i_25', '5v8i_35', '5v8i_37', '4v7x_31', '4v7x_33', '4v7x_46', '4v7x_53', '4v7x_54', '4lf5_5', '4lf5_6', '4ujd_2', '4ujd_3', '4wqr_1', '4wqr_3', '4wqr_9', '4wqr_20', '4wqr_29', '4wqr_34', '4wfn_1', '4wfn_8', '4wfn_10', '4wfn_11', '6gxn_11', '4v5p_11', '4v5p_12', '4v5p_20', '4v5p_25', '4v5p_38', '4v5p_46', '6fyx_2', '5l3p_5', '5l3p_10', '4v6n_4', '5el6_7', '5el6_12', '5el6_22', '5el6_33', '5el6_39', '5el6_56', '6n9f_30', '6n9f_35', '5nd8_6', '5nd8_9', '5nd8_10', '6c4i_17', '5o60_9', '1p9x_4', '4v4r_4', '4v4r_5', '4v4r_7', '4wqf_2', '4wqf_3', '4wqf_6', '4wqf_12', '4wqf_21', '4wqf_54', '4wqf_55', '4dr4_6', '4v7l_12', '4v7l_22', '4v7l_27', '4v7l_41', '2d3o_5', '2d3o_14', '5li0_6', '5li0_7', '5li0_12', '4ybb_22', '4ybb_52', '5e81_8', '5e81_20', '5e81_29', '5e81_32', '5e81_36', '5e81_50', '5kpw_23', '3jaj_6', '3jaj_7', '3jaj_8', '3jaj_9', '3jaj_12', '3j7q_1', '3j7q_2', '5ndj_21', '5ndj_30', '5ndj_37', '5ndj_52', '4yhh_4', '6gz4_4', '6gz4_11', '6gz4_15', '3ccm_14', '5h1s_4', '5lza_10', '6dzi_18', '6dzi_19', '6gqb_1', '1qvf_20', '5apo_5', '5ju8_2', '2j28_2', '1q7y_21', '4v9a_4', '4v9a_36', '4wt1_8', '4wt1_20', '4wt1_29', '4wt1_31', '4wt1_34', '4wt1_47', '3jcn_19', '4v56_37', '4v56_45', '5uym_15', '6ftg_7', '6ftg_9', '4v8o_8', '4v8o_12', '6n1d_16', '6n1d_17', '6d9j_4', '6d9j_9', '6d9j_12', '6gsj_6', '6gsj_30', '6gsj_47', '6gsj_51', '1vqo_22', '4tub_7', '4tub_13', '4tub_21', '3j3v_1', '3j3v_2', '4v8c_11', '4v8c_23', '4v8c_38', '5afi_14', '1i94_8', '5wdt_16', '6gqv_2', '6gqv_3', '4ji2_6', '4lsk_2', '4lsk_35', '4lsk_48', '4lsk_54', '3ccu_14', '4v9m_6', '4v9m_28', '4v9m_36', '5ady_1', '4v8u_1', '4v8u_10', '4v8u_33', '4v8u_37', '4v8u_48', '4v8u_51', '5j30_3', '5j30_10', '5j30_11', '5j30_13', '3j2b_6', '6hcj_2', '1yj9_15', '1yj9_20', '1s72_21', '1ibl_13', '5lzc_8', '6hcf_5', '6dzk_3', '2qex_3', '4b3s_10', '3jah_1', '3jah_5', '3jah_12', '5dox_28', '5dox_29', '5dox_33', '5dox_36', '4v70_3', '2vqf_1', '3jbv_2', '3jbv_13', '3jbv_21', '3jbv_25', '5i4l_4', '5i4l_6', '4nxn_6', '4u1v_6', '4u1v_9', '5fcj_1', '5fcj_4', '5fcj_9', '5fcj_23', '6fti_7', '6fti_9', '4wu1_8', '4wu1_29', '4wu1_40', '4wu1_43', '4wu1_47', '4wu1_49', '4v8a_3', '4v8a_12', '4v8a_24', '4v8a_26', '5czp_4', '5czp_10', '5czp_11', '5czp_13', '3pio_1', '3pio_6', '1i96_8', '4v8m_3', '4v8m_7', '4fb0_2', '5aj0_1', '5aj0_14', '1vqm_20', '6hrm_6', '5u9f_11', '5lzw_11', '4v54_37', '4v54_43', '3g71_13', '3g71_16', '4v9c_6', '4v9c_8', '4v9c_12', '4v9c_42', '4v9c_44', '4v9c_48', '4v9c_51', '1q82_20', '4w29_1', '4w29_24', '4w29_31', '4lf7_8', '3j79_9', '3j79_10', '4v7z_27', '4v7z_29', '4v7z_42', '4v7z_46', '3oto_4', '3oto_5', '3oto_7', '4v83_1', '4v83_46', '4w2g_25', '4w2g_31', '4w2g_34', '4w2g_48', '4v6t_5', '4v6t_15', '4v6t_31', '4v5j_1', '4v5j_30', '4v5j_33', '4v5j_48', '4v5j_50', '4v5j_53', '4y4o_39', '4y4o_40', '4y4o_43', '4y4o_45', '4y4o_59', '4v91_3', '4u26_5', '4u26_7', '4u26_11', '4u26_34', '4u26_50', '4u26_51', '4v5f_1', '4v5f_32', '4v5f_44', '4v5f_47', '4yy3_1', '4v7v_3', '4v7v_7', '4v7v_19', '4v7v_52', '5jup_9', '4v4h_17', '4v4h_32', '4v4h_40', '4v4h_47', '5dm6_2', '5dm6_4', '5dm6_13', '4wce_3', '5vp2_35', '5vp2_38', '5vp2_41', '5vp2_52', '5vp2_55', '5x8t_13', '4v4p_5', '5j88_20', '5j88_46', '6mte_1', '4e8k_2', '5el4_7', '5el4_29', '5el4_42', '5el4_46', '5el4_51', '1n34_3', '1n34_4', '1n34_5', '1n34_7', '3gx7_1', '4v7b_7', '4v7b_12', '4v7b_17', '4v7b_33', '5j4c_12', '5j4c_39', '5j4c_51', '5j4c_52', '5it7_1', '5mei_5', '5mei_7', '4v5r_12', '4v5r_13', '4v5r_23', '4v5r_32', '4v5r_39', '6n8j_4', '6n8j_7', '6n8j_8', '4v6l_6', '4v6l_12', '4v6l_14']

'''11/09/19 liste des A-minor entre deux chaines '''
deux_chaines = [('4wqy', 9), ('4wqy', 29), ('5mdv', 30), ('5zeu', 3), ('5wit', 2), ('5wit', 28), ('4v5c', 17), ('4v5c', 49), ('3j9w', 13), ('5e7k', 12), ('5e7k', 25), ('6b4v', 2), ('5ib8', 14), ('5ib8', 28), ('5mdz', 32), ('6bz8', 17), ('6bz8', 29), ('4wqu', 15), ('4wqu', 39), ('5dfe', 16), ('6cfj', 1), ('6cfj', 27), ('5mgp', 13), ('4v7k', 25), ('4v7k', 27), ('1xmo', 4), ('4wpo', 20), ('4wpo', 23), ('4wpo', 43), ('5f8k', 15), ('4v6e', 15), ('4v6e', 19), ('4v6e', 23), ('2gcs', 1), ('4woi', 53), ('4woi', 62), ('4v9r', 6), ('4v9r', 33), ('4v67', 51), ('4v67', 60), ('5hcp', 5), ('5hcp', 17), ('5lmv', 5), ('4v79', 8), ('1vy4', 24), ('1vy4', 41), ('1vy4', 60), ('4v51', 19), ('4v51', 54), ('3j0o', 2), ('4tue', 21), ('4tue', 48), ('4v8d', 46), ('4v8d', 49), ('2ho7', 1), ('5w4k', 2), ('5w4k', 27), ('4v9j', 28), ('3jce', 20), ('4yzv', 24), ('4yzv', 52), ('1vy6', 35), ('4meh', 1), ('4wq1', 28), ('4wq1', 42), ('5hcr', 5), ('5lmt', 7), ('6gaz', 2), ('5lzd', 17), ('4lfz', 18), ('4lfz', 43), ('4x62', 5), ('4wt8', 20), ('4wt8', 50), ('4v9h', 4), ('3g8s', 1), ('3g8s', 2), ('3g8s', 3), ('3g8s', 4), ('6cc3', 1), ('6cc3', 2), ('3l3c', 1), ('3l3c', 2), ('3l3c', 4), ('3l3c', 5), ('4v8f', 44), ('4v8f', 48), ('3oxm', 1), ('3t1y', 2), ('4v8j', 25), ('4v8j', 54), ('4v9d', 14), ('4v9d', 42), ('6boh', 30), ('6boh', 35), ('4lnt', 26), ('4lnt', 60), ('6q9a', 6), ('3j9y', 15), ('1hmh', 1), ('4v5a', 40), ('5mrf', 11), ('6nd6', 2), ('6nd6', 27), ('3owz', 1), ('3owz', 2), ('6bz6', 16), ('6bz6', 26), ('4v6g', 12), ('4v6g', 51), ('6cae', 1), ('6cae', 2), ('6cae', 26), ('6cae', 32), ('6mtb', 5), ('5j4d', 33), ('5j4d', 36), ('4y4p', 2), ('4y4p', 31), ('6gsk', 9), ('6gsk', 20), ('4tuc', 20), ('4tuc', 46), ('3oxe', 1), ('4v8n', 17), ('4v8n', 28), ('4v8n', 57), ('3deg', 2), ('5uyl', 17), ('5h5u', 32), ('6hcq', 3), ('5lzt', 7), ('5lzx', 6), ('4v9l', 15), ('4v9l', 28), ('6gwt', 7), ('4v8b', 46), ('4v8b', 48), ('3b4a', 1), ('6enj', 8), ('2h0x', 1), ('5vpo', 30), ('5vpo', 45), ('5ndk', 17), ('5ndk', 20), ('4x66', 3), ('6bmd', 1), ('6bmd', 2), ('3j0q', 2), ('4lel', 22), ('4lel', 45), ('4wzo', 10), ('4wzo', 23), ('4wzo', 24), ('3jag', 11), ('4v61', 8), ('5d8b', 7), ('6enf', 13), ('5o61', 26), ('4v7a', 15), ('5el7', 15), ('5el7', 28), ('4v5q', 41), ('4v5q', 53), ('5kpv', 2), ('5jte', 12), ('4kzz', 2), ('4v7m', 32), ('4v7m', 36), ('4v7m', 39), ('4v7m', 42), ('4v4s', 8), ('4ypb', 21), ('4ypb', 49), ('1xmq', 4), ('4v5e', 40), ('4v5e', 54), ('4w2h', 36), ('4w2h', 39), ('5hd1', 6), ('5hd1', 18), ('4wra', 36), ('4wra', 44), ('5t2a', 5), ('4w4g', 22), ('4w4g', 49), ('6buw', 14), ('6buw', 26), ('5j3c', 16), ('5j3c', 29), ('4v6a', 16), ('5kpx', 2), ('4v5s', 46), ('4v5s', 54), ('4v5s', 61), ('5j4b', 6), ('5j4b', 18), ('6n9e', 3), ('6n9e', 16), ('5el5', 15), ('4w2f', 40), ('4w2f', 55), ('6gxm', 9), ('3j78', 5), ('5tcu', 4), ('4wro', 38), ('4v4i', 6), ('4tzz', 1), ('4tzz', 2), ('6ha1', 11), ('4v90', 5), ('4v5g', 25), ('4v5g', 37), ('4v5g', 43), ('4v5g', 47), ('4wsm', 34), ('3jaq', 2), ('3b4c', 1), ('3jbo', 10), ('6q95', 5), ('4v9n', 25), ('4v9n', 44), ('5lzz', 8), ('5wfk', 14), ('2uxd', 2), ('5lzv', 6), ('4v9b', 23), ('4v9b', 26), ('4tua', 22), ('4tua', 49), ('5u9g', 18), ('5ibb', 16), ('5ibb', 28), ('5lmr', 6), ('5wfs', 12), ('3j2c', 1), ('3j2c', 3), ('4v63', 33), ('4v63', 42), ('2il9', 1), ('2il9', 2), ('4v8x', 15), ('4v8x', 41), ('6o97', 1), ('6o97', 3), ('6o97', 29), ('3g9c', 1), ('3g9c', 2), ('3g9c', 3), ('3g9c', 5), ('3jai', 10), ('5doy', 1), ('5doy', 26), ('5o2r', 14), ('4tud', 21), ('4tud', 48), ('3oxb', 1), ('6gsl', 4), ('6gsl', 29), ('6gsl', 30), ('6gsl', 42), ('6bok', 20), ('6bok', 29), ('5lzs', 8), ('4v50', 13), ('4v50', 38), ('6enu', 15), ('4p70', 23), ('4p70', 45), ('4zer', 13), ('4zer', 45), ('4v9k', 16), ('3jcd', 3), ('6ndk', 12), ('5hau', 12), ('5hau', 57), ('2ho6', 1), ('4v8e', 30), ('4v8e', 48), ('4gkj', 2), ('4v9s', 6), ('4v9s', 30), ('4meg', 1), ('4wzd', 15), ('4wzd', 27), ('5j8b', 24), ('4v8q', 22), ('1vy5', 21), ('1vy5', 38), ('1vy5', 54), ('1vy5', 55), ('4v78', 11), ('5hcq', 5), ('5hcq', 17), ('5uq7', 4), ('5kcr', 17), ('5we6', 15), ('6of1', 1), ('6of1', 2), ('6of1', 30), ('4v6d', 14), ('4v6d', 22), ('4v7j', 30), ('4v7j', 32), ('5x8p', 3), ('5x8p', 30), ('4v4t', 11), ('4l71', 21), ('4l71', 46), ('5mre', 11), ('4v95', 24), ('4v95', 54), ('5u4i', 9), ('5mdw', 31), ('5jbh', 3), ('6nd5', 3), ('6nd5', 30), ('4l47', 23), ('4l47', 52), ('3ox0', 1), ('3ox0', 2), ('4v87', 49), ('6ha8', 11), ('5aa0', 5), ('5wf0', 14), ('3j9z', 14), ('4wsd', 27), ('4wsd', 39), ('5kps', 7), ('4v6f', 35), ('4v6f', 42), ('1xnr', 4), ('5we4', 16), ('3g96', 1), ('3g96', 2), ('3g96', 4), ('3g96', 5), ('4v4z', 11), ('5zeb', 4), ('5wnv', 4), ('4v5l', 29), ('4z8c', 5), ('4z8c', 19), ('3oww', 1), ('5mdy', 29), ('2oiu', 1), ('2oiu', 2), ('5ib7', 18), ('5ib7', 31), ('5ib7', 32), ('6bz7', 14), ('6bz7', 23), ('4v7p', 14), ('4v7p', 28), ('4v97', 23), ('4v97', 51), ('5nwy', 12), ('1vvj', 21), ('1vvj', 46), ('4v9i', 18), ('4v9i', 54), ('5mmm', 23), ('3jcj', 19), ('4wr6', 37), ('4zsn', 22), ('4zsn', 49), ('5mlc', 5), ('3j0l', 2), ('5me1', 5), ('5lmu', 4), ('6gaw', 7), ('5iqr', 5), ('5iqr', 14), ('1vy7', 36), ('4v9q', 40), ('4v9q', 42), ('4z3s', 3), ('4z3s', 33), ('5lze', 13), ('4v68', 23), ('4p6f', 20), ('4p6f', 45), ('4v4j', 9), ('4w2i', 41), ('4w2i', 56), ('4lt8', 24), ('4lt8', 54), ('5wis', 4), ('5wis', 15), ('4v5d', 8), ('4v5d', 39), ('4v5d', 57), ('5mrc', 11), ('4w2e', 16), ('4wqr', 39), ('6gxn', 9), ('4v5p', 29), ('4v5p', 53), ('5l3p', 9), ('5el6', 31), ('5el6', 45), ('6n9f', 3), ('6n9f', 15), ('6c4i', 14), ('4v4r', 9), ('4wqf', 14), ('4v7l', 5), ('4v7l', 28), ('4v7l', 31), ('5e81', 27), ('5e81', 42), ('5kpw', 3), ('6i7o', 6), ('6i7o', 24), ('5ndj', 13), ('5ndj', 16), ('3t1h', 3), ('5lza', 11), ('5lmq', 7), ('5ju8', 11), ('3j0p', 2), ('2nz4', 1), ('2nz4', 2), ('2nz4', 4), ('2nz4', 5), ('4v9a', 23), ('4v9a', 27), ('5mmi', 2), ('5mmi', 20), ('4wt1', 37), ('5lzu', 6), ('3jcn', 19), ('5uym', 25), ('4v8o', 18), ('6n1d', 13), ('4plx', 1), ('4plx', 4), ('6gsj', 14), ('6gsj', 27), ('4tub', 19), ('4tub', 43), ('3oxd', 1), ('3oxd', 2), ('4v8c', 44), ('5wdt', 15), ('4lsk', 25), ('4lsk', 53), ('5lzy', 3), ('4v9m', 16), ('4v9m', 30), ('4v8u', 42), ('4v8u', 58), ('2h0w', 1), ('6hcj', 3), ('2gcv', 1), ('5lzc', 9), ('3jah', 10), ('3g8t', 1), ('3g8t', 2), ('3g8t', 3), ('3jbn', 10), ('4wu1', 15), ('4wu1', 26), ('3b4b', 1), ('5czp', 15), ('5czp', 26), ('3oxj', 1), ('4v8m', 3), ('4v8m', 5), ('5u9f', 16), ('5lzw', 6), ('4v9c', 21), ('4v9c', 45), ('4w29', 22), ('4v83', 29), ('4v83', 38), ('4w2g', 35), ('4w2g', 38), ('4v5j', 42), ('4v5f', 23), ('4v5f', 48), ('4dr6', 5), ('5vp2', 5), ('5vp2', 16), ('5x8t', 3), ('5x8t', 22), ('4v4p', 1), ('5el4', 14), ('5el4', 25), ('5el4', 26), ('3owi', 1), ('4v7b', 16), ('5j4c', 1), ('5j4c', 27), ('4v5r', 30), ('4v5r', 52)]


