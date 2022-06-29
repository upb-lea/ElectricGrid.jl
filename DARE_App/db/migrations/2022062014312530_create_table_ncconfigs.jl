module CreateTableNcconfigs

import SearchLight.Migrations: create_table, column, columns, pk, add_index, drop_table, add_indices

function up()
  create_table(:ncconfigs) do
    [
      pk()
      column(:name, :string, limit = 100)
      column(:cm, :string, limit = 1000000)
      column(:parameters, :string, limit = 1000000)
    ]
  end

  add_index(:ncconfigs, :name)
end

function down()
  drop_table(:ncconfigs)
end

end
